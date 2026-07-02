#include "source_base.h"

namespace twinkle
{

using f_t = source_base_t::f_t;
using c_t = source_base_t::c_t;

__device__ int source_base_t::prev_src_idx_g
( const int idx_here, const int srcidx ) const
{
	bool skip;
	int idx_prev;
	idx_prev = pool_margin[ srcidx * n_point_max + idx_here].prev_src_idx;
	skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
	while(skip)
	{
		idx_prev = pool_margin[ srcidx * n_point_max + idx_prev ].prev_src_idx;
		skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
	}
	return idx_prev;
} 

__device__ int source_base_t::next_src_idx_g
( const int idx_here, const int srcidx ) const
{
	bool skip;
	int idx_next;
	idx_next = pool_margin[ srcidx * n_point_max + idx_here ].next_src_idx;
	skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
	while(skip)
	{
		idx_next = pool_margin[ srcidx * n_point_max + idx_next ].next_src_idx;
		skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
	}
	return idx_next;
}


__device__ int source_base_t::prev_src_idx_local
( const int idx_here, local_info_t<f_t>& local_info, bool * changed ) const
{
	bool skip;
	int idx_prev;
    int batchidx = local_info.batchidx;

    idx_prev = local_info.new_pts[idx_here - batchidx*n_th].prev_src_idx;
    if(idx_prev >= batchidx*n_th)       // 如果是新点，验证是否skip
    {
        skip = local_info.new_pts[idx_prev - batchidx*n_th].skip;
        if(skip){*changed = true;}
        while(skip)
        {
            idx_prev = local_info.new_pts[idx_prev - batchidx*n_th].prev_src_idx;
            if(idx_prev < batchidx*n_th)
            {
                break;
            }
            skip = local_info.new_pts[idx_prev - batchidx*n_th].skip;
        }
    }
    // else{};  // 如果是旧点，一定不会skip，因为上一个循环改好了
    return idx_prev;
}

__device__ int source_base_t::next_src_idx_local
( const int idx_here, local_info_t<f_t>& local_info, bool * changed ) const
{
	bool skip;
	int idx_next;
    int batchidx = local_info.batchidx;

    idx_next = local_info.new_pts[idx_here - batchidx*n_th].next_src_idx;
    if(idx_next >= batchidx*n_th)
    {
        skip = local_info.new_pts[idx_next - batchidx*n_th].skip;
        if(skip){*changed = true;}
        while(skip)
        {
            idx_next = local_info.new_pts[idx_next - batchidx*n_th].next_src_idx;
            if(idx_next < batchidx*n_th)
            {
                break;          // 如果找到旧点了，它一定没有skip，就是它
            }
            skip = local_info.new_pts[idx_next - batchidx*n_th].skip;
        }
    } // 逻辑见 prev 版本
    return idx_next;
}

__device__ int source_base_t::prev_src_idx_local_g
( const int idx_here, local_info_t<f_t>& local_info, const int srcidx ) const
{
	bool skip;
	int idx_prev;
    int batchidx = local_info.batchidx;
    bool global = false;

    idx_prev = local_info.new_pts[idx_here - batchidx*n_th].prev_src_idx;
    if(idx_prev < batchidx*n_th)
    {    
        skip = local_info.new_pts[idx_here - batchidx*n_th].skip;
        while(skip)
        {
            idx_prev = local_info.new_pts[idx_prev - batchidx*n_th].prev_src_idx;
            if(idx_prev < batchidx*n_th)
            {
                global = true;
                break;
            }
            skip = local_info.new_pts[idx_prev - batchidx*n_th].skip;
        }
    }
    if(global)
    {
        skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
        while(skip)
        {
            idx_prev = pool_margin[ srcidx * n_point_max + idx_prev ].prev_src_idx;
            skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
        }
    }

    return idx_prev;
}

__device__ int source_base_t::next_src_idx_local_g
( const int idx_here, local_info_t<f_t>& local_info, const int srcidx ) const
{
	bool skip;
	int idx_next;
    int batchidx = local_info.batchidx;
    bool global = false;

    idx_next = local_info.new_pts[idx_here - batchidx*n_th].next_src_idx;
    if(idx_next < batchidx*n_th)
    {
        skip = local_info.new_pts[idx_here - batchidx*n_th].skip;
        while(skip)
        {
            idx_next = local_info.new_pts[idx_next - batchidx*n_th].next_src_idx;
            if(idx_next < batchidx*n_th)
            {
                global = true;
                break;
            }
            skip = local_info.new_pts[idx_next - batchidx*n_th].skip;
        }        
    }
    if(global)
    {
        skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
        while(skip)
        {
            idx_next = pool_margin[ srcidx * n_point_max + idx_next ].next_src_idx;
            skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
        }        
    }
    
    return idx_next;
}



__device__ void source_base_t::connect_order
( const c_t* posi_here, const c_t* posi_other, int len_here, int len_other, int* order) const
{
	f_t Distance[6];
	int jbest=0;
	f_t temp;
	c_t temp_c;
	if((len_here==2))
	{
		if(len_other==2)
		{
			// f_t Distance[2];	// 01,10
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );

			if(Distance[0] < Distance[1])
			{
				order[0] = 0;
				order[1] = 1;
			}
			else
			{
				order[0] = 1;
				order[1] = 0;
			}			
		}
		else	// len_other==3
		{
			// f_t Distance[6];	// 01,02,12,10,20,21
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  );
			Distance[2] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  );
			Distance[3] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );
			Distance[4] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );
			Distance[5] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );

			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1.3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
		}
	}
	else		// len_here==3
	{
		if(len_other==2)
		{
			// f_t Distance[6];	// 01X，0X1, 1X0, 10X, X01, X10
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[2] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );
			Distance[3] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );
			Distance[4] = (posi_here[1]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[5] = (posi_here[1]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );

			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1.3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
			order[2] = 2 - (jbest%3);	// 2,1,0,2,1,0

			for(int j=0;j<3;j++)
			{
				if(order[j]==2){order[j]=-1;}	// 2为虚解，标记为未连接
			}
		}
		else	// len_other==3
		{
			// f_t Distance[6];	// 012, 021, 120, 102, 201, 210
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[2]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[2] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );
			Distance[3] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[2]).abs(  );
			Distance[4] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[5] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );

			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1,3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
			order[2] = 2 - (jbest%3);	// 2,1,0,2,1,0
		}
	}
}

__device__ int source_base_t::positive_connect_j34     // return j_positive_candidate
( int* next_idx_out, int* next_j_out, const int Nphys_here, const int Nphys_next,
    const c_t * posi_here_positive, const c_t * posi_next, const int here_idx, const int next_idx) const
{
    int j_positive_candidate = -1;
    if(Nphys_here==3)
    {
        next_idx_out[3] = -13;
        next_j_out[3] = -13;  
    }
    if(Nphys_here==Nphys_next)
    {
        if(Nphys_here==3)
        {
            // output
            next_idx_out[4] = next_idx;
            next_j_out[4] = 4;
        }
        else    // Nphys_here == 5, 两个比较大小
        {
            int order_e[2];
            connect_order(posi_here_positive,posi_next,2,2,order_e);
            // output
            // j_here_positive = j_next = {4,3}
            next_idx_out[4] = next_idx;
            next_idx_out[3] = next_idx;
            next_j_out[4] = (order_e[0]==0 ? 4 : 3);
            next_j_out[3] = (order_e[1]==1 ? 3 : 4);
        }
    }
    else    // Nphys_here!=Nphys_next
    {
        if(Nphys_here > Nphys_next)     // here: 2, next: 1
        {
            if((posi_here_positive[0]-posi_next[0]).norm2(  ) < (posi_here_positive[1]-posi_next[0]).norm2(  )) // 连近的，现在00是近的
            {
                // output
                next_idx_out[4] = next_idx;
                next_j_out[4] = 4;
                next_idx_out[3] = here_idx;        // connect to here
                j_positive_candidate = 3;       // j3 还没连
            }
            else    // 10是近的
            {
                // output
                next_idx_out[3] = next_idx;
                next_j_out[3] = 4;
                next_idx_out[4] = here_idx;        // connect to here                            
                j_positive_candidate = 4;       // j4 还没连
            }
        }
        else            // here: 1, next: 2
        {
            if((posi_here_positive[0]-posi_next[0]).norm2(  ) < (posi_here_positive[0]-posi_next[1]).norm2(  )) // 连近的，现在00是近的
            {
                // output
                // pt_here.next_idx[4] = next_idx;
                // pt_here.next_j[4] = 4;
                next_idx_out[4] = next_idx;
                next_j_out[4] = 4;
            }            
            else    // 01是近的
            {
                // output
                // pt_here.next_idx[4] = next_idx;
                // pt_here.next_j[4] = 3;   
                next_idx_out[4] = next_idx;
                next_j_out[4] = 3;   
            }            
        }
    }
    return j_positive_candidate;
}

__device__ int source_base_t::negative_connect_j012    // return j_negative_candidate
( int* next_idx_out, int* next_j_out, const int Nphys_here, const int Nphys_prev,
    const c_t * posi_here_negative, const c_t * posi_prev, const int here_idx, const int prev_idx) const
{
    int order_o[3];
    int j_negative_candidate = -1;

    if(Nphys_here==3)
    {
        next_idx_out[2] = -14;
        next_j_out[2] = -14;  
    }
    if(Nphys_here==Nphys_prev)
    {
        if(Nphys_here==3)
        {
            connect_order(posi_here_negative, posi_prev,2,2,order_o);
            // output
            for(int ii=0;ii<2;ii++)
            {
                next_idx_out[ii] = prev_idx;
                next_j_out[ii] = order_o[ii];
            }           
        }
        else
        {
            connect_order(posi_here_negative, posi_prev,3,3,order_o);
            for(int ii=0;ii<3;ii++)
            {
                next_idx_out[ii] = prev_idx;
                next_j_out[ii] = order_o[ii];
            }                        
        }
    }
    else
    {
        if(Nphys_here > Nphys_prev) // here: 3, prev: 2 
        {
            connect_order(posi_here_negative,posi_prev,3,2,order_o);      // order is 0,1 or -1 (cross candidate)
            // output
            for(int ii=0;ii<3;ii++)
            {
                if(order_o[ii]>=0)  // real
                {
                    next_idx_out[ii] = prev_idx;
                    next_j_out[ii] = order_o[ii];                                
                }
                else    // cross candidate
                {
                    next_idx_out[ii] = here_idx;    // connect to here
                    j_negative_candidate = ii;
                }
            }                          
        }
        else        // here: 2, prev: 3
        {

            connect_order(posi_here_negative,posi_prev,2,3,order_o);

            // output
            for(int ii=0;ii<2;ii++)
            {
                next_idx_out[ii] = prev_idx;
                next_j_out[ii] = order_o[ii];
            }                        
        }
    }
    return j_negative_candidate;
}

__device__ bool source_base_t::cross_j34              // 可能修改 012
( int* next_j_out, const int j_positive_candidate,
const c_t * posi_here_negative, const c_t * posi_next_negative, const int here_idx) const
{
    // int temp1;
    int order_c[3];     // c means cross caustic
    bool fail = false;

    connect_order(posi_here_negative,posi_next_negative,3,2,order_c);
    for(int ii=0;ii<3;ii++)
    {
        if(order_c[ii]==-1)     // 连接会有一个落单的，即离每个nextnegative都很远，这个点被positive cross连接
        {
            next_j_out[j_positive_candidate] = ii;
        }
    }
    return fail;
}

__device__ bool source_base_t::cross_j012               // 可能修改34
( int* next_j_out, const int j_negative_candidate,
const c_t * posi_here_positive, const c_t * posi_prev_positive, const int here_idx) const
{
    // int temp0;
    bool fail = false;

    // posi_here_positive[2] connect to posi_prev_positive[1]
    if((posi_here_positive[0]-posi_prev_positive[0]).norm2(  ) < (posi_here_positive[1]-posi_prev_positive[0]).norm2(  ))
    {
        next_j_out[j_negative_candidate] = 3;
    }
    else    // posi_here_positive[0]落单
    {
        next_j_out[j_negative_candidate] = 4;
    }
    return fail;
}


// j_here_negative = j_prev = {0,1,2} or {0,1,-1}
__device__ bool source_base_t::connect_prev_j012
( int* next_idx_out, int* next_j_out,
const img_t * imgs_here, const img_t * imgs_prev, const int here_idx, const int prev_idx, 
const int Nphys_here, const int Nphys_prev) const
{
    c_t posi_here_negative[3];
    c_t posi_prev[3];
    int j_negative_candidate;
    bool fail = false;

    posi_here_negative[0] = imgs_here[0].position;
    posi_here_negative[1] = imgs_here[1].position;
    if(Nphys_here == 5)
    {
        posi_here_negative[2] = imgs_here[2].position;
    }
    else
    {
        posi_here_negative[2] = 0;             
    }
    posi_prev[0] = imgs_prev[0].position;
    posi_prev[1] = imgs_prev[1].position;
    if(Nphys_prev == 5)
    {
        posi_prev[2] = imgs_prev[2].position;
    }
    else
    {
        posi_prev[2] =0;
    }

    j_negative_candidate = negative_connect_j012
    ( next_idx_out, next_j_out, Nphys_here, Nphys_prev,
    posi_here_negative, posi_prev, here_idx, prev_idx);      // only change next_idx_out[0,1,2], next_j_out[0,1,2]

    if(j_negative_candidate>=0 && j_negative_candidate<order-1)
    {
        c_t posi_prev_positive[1];   // Nphys==3
        c_t posi_here_positive[2];
        posi_prev_positive[0] = imgs_prev[4].position;
        posi_here_positive[0] = imgs_here[4].position;
        posi_here_positive[1] = imgs_here[3].position;
        fail = cross_j012
        ( next_j_out, j_negative_candidate,
        posi_here_positive, posi_prev_positive, here_idx);
    }
    return fail;
}

// j_here_positive = j_next = {4,3} or {4, -1}
__device__ bool source_base_t::connect_next_j34
( int* next_idx_out, int* next_j_out,
const img_t * imgs_here, const img_t * imgs_next, const int here_idx, const int next_idx,
const int Nphys_here, const int Nphys_next) const
{
    c_t posi_here_positive[2];
    c_t posi_next[2];
    int j_positive_candidate;
    bool fail = false;

    posi_here_positive[0] = imgs_here[4].position;
    if(Nphys_here == 5)
    {
        posi_here_positive[1] = imgs_here[3].position;
    }
    else
    {
        posi_here_positive[1] = 0;      // 保险起见的初始化
    }
    posi_next[0] = imgs_next[4].position;
    if(Nphys_next == 5)
    {
        posi_next[1] = imgs_next[3].position;
    }
    else
    {
        posi_next[1] = 0;
    }            

    // positive parity, 1c1 or 2c2
    j_positive_candidate = positive_connect_j34
    ( next_idx_out, next_j_out, Nphys_here, Nphys_next, 
    posi_here_positive, posi_next, here_idx, next_idx);      // only change next_idx_out[3,4], next_j_out[3,4]
    // 检查跨焦散
    // positive parity connect to here
    if(j_positive_candidate>=0 && j_positive_candidate<order-1 )
    {
        c_t posi_next_negative[2];    // Nphys==3
        c_t posi_here_negative[3];
        posi_next_negative[0] = imgs_next[0].position;
        posi_next_negative[1] = imgs_next[1].position;
        posi_here_negative[0] = imgs_here[0].position;
        posi_here_negative[1] = imgs_here[1].position;
        posi_here_negative[2] = imgs_here[2].position;
        fail = cross_j34
        ( next_j_out, j_positive_candidate,
        posi_here_negative, posi_next_negative, here_idx );
    }
    return fail;
}

}; // namespace twinkle
