# Twinkle
A GPU-based High-efficiency binary microlensing code.

If you use Twinkle (or part of Twinkle) in your research work, we request citing our paper: Wang, S. et al., 2024, in prep.

If you incorporate Twinkle  (or part of Twinkle) in your code, we request specifying it in the README file (or equivalent), including the GitHub link ( https://github.com/AsterLight0626/Twinkle ) and asking the user to cite our paper: Wang, S. et al., 2024, in prep.

# 使用教程 
推荐使用CUDA进行计算。使用提供的 MakeFile 进行编译

## 批量计算放大率
使用 /src/demo.cpp 做演示，提供相关参数即可输出放大率结果

## 拟合观测的光变曲线
使用 /src/MCMCdemo.cpp 做演示，此处的 MCMC 算法建议替换为其他更成熟的算法

## 计算临界曲线和焦散线
使用 /src/caus_demo.cpp 做演示，输出 causcrit.txt

## 使用提供的 jupyter notebook 文件进行结果可视化
