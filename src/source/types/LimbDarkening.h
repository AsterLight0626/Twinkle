#pragma once

#include <array>
#include <cmath>

namespace twinkle
{
////////////////////////////////////////////////////////////
// 自适应辛普森积分误差单元结构体
// 所有成员均为 public 以便调试
template< class f_T >
struct ErrorUnit_t
{
    // 核心数据数组
    std::array<f_T, 5> phi;     // [L, LC, C, RC, R]
    std::array<f_T, 5> val;     // [val_L, val_LC, val_C, val_RC, val_R]
    std::array<int, 5> Ncross;  // [nc_L, nc_LC, nc_C, nc_RC, nc_R]
    std::array<f_T, 5> raw_Mag; // 2603 新增：用于 monotest

    // 访问函数（替代引用成员）
    f_T& L() { return phi[0]; }
    const f_T& L() const { return phi[0]; }
    f_T& LC() { return phi[1]; }
    const f_T& LC() const { return phi[1]; }
    f_T& C() { return phi[2]; }
    const f_T& C() const { return phi[2]; }
    f_T& RC() { return phi[3]; }
    const f_T& RC() const { return phi[3]; }
    f_T& R() { return phi[4]; }
    const f_T& R() const { return phi[4]; }

    f_T& val_L() { return val[0]; }
    const f_T& val_L() const { return val[0]; }
    f_T& val_LC() { return val[1]; }
    const f_T& val_LC() const { return val[1]; }
    f_T& val_C() { return val[2]; }
    const f_T& val_C() const { return val[2]; }
    f_T& val_RC() { return val[3]; }
    const f_T& val_RC() const { return val[3]; }
    f_T& val_R() { return val[4]; }
    const f_T& val_R() const { return val[4]; }

    // 其他成员
    f_T h;                       // 区间长度 R-L
    f_T SimpsonC, SimpsonL, SimpsonR;  // 辛普森积分值
    f_T Error;                   // 误差估计
    int depth;                   // 递归深度
    f_T ub, lb;                  // 上下界
    f_T* phi_hidden = nullptr;             // hidden error

    static constexpr double coeff_ec = 0.1;      // 2603 新增：cross error

    // 构造函数
    ErrorUnit_t(f_T left, f_T right) 
        : phi{}
        , val{}
        , Ncross{}
        , h(right - left)
        , SimpsonC(0)
        , SimpsonL(0)
        , SimpsonR(0)
        , Error(0)
        , depth(0)
        , ub(0)
        , lb(0)
    {
        phi[0] = left;
        phi[4] = right;
        for(int i = 1; i < 4; ++i) {
            phi[i] = left + i * h / 4.0;
        }
        val.fill(static_cast<f_T>(0));
        Ncross.fill(-1);
        raw_Mag.fill(-1.);
    }

    ErrorUnit_t() 
        : ErrorUnit_t(0.0, 1.0)  // 委托构造
    {}

    // 比较运算符（用于优先队列等容器）
    bool operator<(const ErrorUnit_t& other) const
    {
        // return Error/sqrt(SimpsonL+SimpsonR) < other.Error/sqrt(other.SimpsonL+other.SimpsonR);
        return Error < other.Error;
    }
    
    
    // 计算左右子区间的辛普森积分
    void SimpsonLR()
    {
        SimpsonL = h / 12 * (val_L() + 4 * val_LC() + val_C());
        SimpsonR = h / 12 * (val_C() + 4 * val_RC() + val_R());
    }
    
    
    // 分裂为左右两个子单元
    std::pair<ErrorUnit_t, ErrorUnit_t> Split()
    {
        // 左子单元: [L, C]
        ErrorUnit_t euL(phi[0], phi[2]);
        euL.val[0] = val[0];  // val_L
        euL.val[4] = val[2];  // val_R = val_C
        euL.val[2] = val[1];  // val_C = val_LC
        euL.Ncross[0] = Ncross[0];
        euL.Ncross[4] = Ncross[2];
        euL.Ncross[2] = Ncross[1];
        euL.SimpsonC = SimpsonL;
        euL.depth = depth + 1;
        euL.ub = ub;  // 继承容差，后续会重新分配
        euL.phi_hidden = phi_hidden;
        
        // 右子单元: [C, R]
        ErrorUnit_t euR(phi[2], phi[4]);
        euR.val[4] = val[4];  // val_R
        euR.val[0] = val[2];  // val_L = val_C
        euR.val[2] = val[3];  // val_C = val_RC
        euR.Ncross[4] = Ncross[4];
        euR.Ncross[0] = Ncross[2];
        euR.Ncross[2] = Ncross[3];
        euR.SimpsonC = SimpsonR;
        euR.depth = depth + 1;
        euR.ub = ub;  // 继承容差，后续会重新分配
        euR.phi_hidden = phi_hidden;
        
        return {euL, euR};
    }
    
    // 计算上下界
    void bound()
    {
        ub = (val_L() + val_LC() + val_C() + val_RC()) / 4 * h;
        lb = (val_LC() + val_C() + val_RC() + val_R()) / 4 * h;
    }

    // 计算误差（要求已给定 val_RC, val_LC, SimpsonC）
    void GetError()
    {
        SimpsonLR();
        Error = std::abs(SimpsonL + SimpsonR - SimpsonC);
        bound();
        // 添加额外误差项 cross error
        if(Ncross[0] != Ncross[2] || Ncross[2] != Ncross[4])
        {
            // Error += (ub - lb) * (R() - L());
            Error += coeff_ec * (ub - lb);
        }
        // hidden error
        if(phi_hidden)
        {
            if(phi_hidden[4]>0 && Ncross[0] == 0 && Ncross[4] == 0)
            {
                for(int phi_i=0;phi_i<int(phi_hidden[4]);phi_i++)
                {
                    if(phi[0] < phi_hidden[phi_i] && phi_hidden[phi_i] < phi[4])
                    {
                        // printf("hidden error added, phiL: %.16f, phiR: %.16f, depth: %d\n",phi[0], phi[4],depth);
                        Error += (ub - lb);
                    }
                }
            }
        }
    }

    // 数值稳定性检验（新增，来自 Python monotest）
    void monotest(const ErrorUnit_t& parent, f_T phi_self, f_T val_self_raw, 
                  f_T& val_out, f_T& mag_out)
    {

        f_T mag_interp = (parent.Mag[0] + parent.Mag[4]) / 2.0;
        f_T val_interp = mag_interp * (1 - phi_self * phi_self);
        
        // 如果偏差过大，使用插值替代
        if(val_self_raw > parent.val[0] * 10 || val_self_raw < parent.val[4] * 0.1) {
            val_out = val_interp;
            mag_out = mag_interp;
        } else {
            val_out = val_self_raw;
            mag_out = val_self_raw / (1 - phi_self * phi_self);
        }
    }
};

}  // namespace twinkle