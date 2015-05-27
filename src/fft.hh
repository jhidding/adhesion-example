#pragma once
#include <fftw3.h>
#include <memory>
#include <complex>
#include <vector>

#include "boxparam.hh"
#include <iostream>

template <typename T>
class FFTW_allocator: public std::allocator<T>
{
    public:
        typedef T		    value_type;
        typedef T *	 	    pointer;
        typedef T &	 	    reference;
        typedef T const * 	const_pointer;
        typedef T const &	const_reference;
        typedef size_t 		size_type;
        typedef ptrdiff_t 	difference_type;

        pointer allocate(size_t n, std::allocator<void>::const_pointer hint = 0)
        {
            if (hint != 0)
                fftw_free(hint);

            return reinterpret_cast<T *>(fftw_malloc(n * sizeof(T)));
        }

        void deallocate(pointer p, size_t n)
        {
            fftw_free(p);
        }
};

class FFT3
{
public:
    using c64 = std::complex<double>;
    std::vector<c64, FFTW_allocator<c64>>
                        fourier_space, real_space;

private:
    BoxParam            box;
    fftw_plan	        d_plan_fwd, d_plan_bwd;

public:
    FFT3(BoxParam const &box_):
        box(box_),
        fourier_space(box_.size),
        real_space(box_.size)
    {
        int N = static_cast<int>(box.N);
        int shape[3] = {N, N, N};

    	d_plan_fwd = fftw_plan_dft(3, shape,
    		reinterpret_cast<fftw_complex *>(real_space.data()),
    		reinterpret_cast<fftw_complex *>(fourier_space.data()),
            FFTW_FORWARD, FFTW_ESTIMATE);

    	d_plan_bwd = fftw_plan_dft(3, shape,
    		reinterpret_cast<fftw_complex *>(fourier_space.data()),
    		reinterpret_cast<fftw_complex *>(real_space.data()),
            FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    void forward_transform()
    {
        fftw_execute(d_plan_fwd);
    }

    void backward_transform()
    {
        fftw_execute(d_plan_bwd);
        for (c64 &z : real_space) z /= box.size;
    }
};
