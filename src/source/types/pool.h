
template< class dat_T >
struct    pool_t
{
    ////////// Data //////////
    using dat_t = dat_T;
    int     stride;
    int      n_pts;
    dat_T *  dat_d;
    dat_T *  dat_h;

    ////////// Host-side interfaces //////////
    __host__ pool_t(  ) : dat_d( nullptr ) ,
                          dat_h( nullptr ) {  };
    
    __host__ void setup_mem
    ( device_t  & f_dev, const int & size,
      const int & stride = 1 )
    {
        this  ->stride = stride;
        this  ->n_pts =    size;
        dat_h = f_dev.malloc_host  < dat_T > ( size );
        dat_d = f_dev.malloc_device< dat_T > ( size );
        return;
    };
    __host__ void free( device_t & f_dev )
    {
        if( dat_d != nullptr )
            f_dev.free_device( dat_d );
        if( dat_h != nullptr )        
            f_dev.free_host  ( dat_h );
        return;
    };
    __host__ void cp_h2d( device_t & f_dev )
    {
        return  f_dev.cp( dat_d, dat_h, n_pts );
    };
    __host__ void cp_d2h( device_t & f_dev )
    {
        return  f_dev.cp( dat_h, dat_d, n_pts );
    };
    
    ////////// Device-side interfaces //////////
    __device__ __forceinline__  dat_T * dispatch
    ( int & len, const int & idx = -1 ) const
    {
        len  = stride;
        auto * dat = const_cast < dat_T * > ( dat_d );
        return dat + len * ( idx < 0 ? blockIdx.x : idx );
    };

    __host__ __device__ __forceinline__
    dat_T & operator [  ] ( const int & idx ) const
    {
        return const_cast< dat_T * >( dat_d )[ idx ];
    };
};
