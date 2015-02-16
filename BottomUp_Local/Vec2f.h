//
//  Vec2f.h
//  Iterative_Particles
//
//  Created by Sean R Purser-Haskell on 1/27/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#ifndef Iterative_Particles_Vec2f_h
#define Iterative_Particles_Vec2f_h

#include <cmath>
#include <cstdlib>
#include <cassert>

typedef unsigned int uint32;

#define PI						3.141592

template<typename T>
struct VecBase2
{
    typedef T Type;
    static const unsigned mComponentCount=2;
    
    VecBase2()
    {
        for(unsigned i=0;i<mComponentCount;++i)
            mComponents[i]=0;
    }
    
    T				&GetComponent(unsigned index)
    {
        return mComponents[index];
    }
    
    T		   const&GetComponentConst(unsigned index)const
    {
        return mComponents[index];
    }
    
    union
    {
        struct
        {
            T x,y;
        };
        struct
        {
            T width,height;
        };
        struct
        {
            T row,column;
        };
        struct
        {
            T longitude, latitude;
        };
        struct
        {
            T u,v;
        };
        
        T mComponents[2];
    };
};



template<typename Base>
struct Vec : public Base
{
    typedef typename Base::Type T;
    
    Vec()
    {
    }
    explicit Vec(T x)
    {
        assert(this->mComponentCount==1);
        
        this->mComponents[0]=x;
    }
    Vec(T x,T y)
    {
        assert(this->mComponentCount==2);
        
        this->mComponents[0]=x;
        this->mComponents[1]=y;
    }
    Vec(T x,T y,T z)
    {
        assert(this->mComponentCount==3);
        
        this->mComponents[0]=x;
        this->mComponents[1]=y;
        this->mComponents[2]=z;
    }
    Vec(T x,T y,T z,T w)
    {
        assert(this->mComponentCount==4);
        
        this->mComponents[0]=x;
        this->mComponents[1]=y;
        this->mComponents[2]=z;
        this->mComponents[3]=w;
    }
    
    bool operator==(Vec const&o)const
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
        {
            if(this->mComponents[i]!=o.mComponents[i])
                return false;
        }
        return true;
    }
    
    bool operator!=(Vec const&o)const
    {
        return !((*this)==o);
    }
    
    bool operator<(Vec const&o)const
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
        {
            if(this->mComponents[i]!=o.mComponents[i])
                return this->mComponents[i]<o.mComponents[i];
        }
        return false;
    }
    bool operator>(Vec const&o)const
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
        {
            if(this->mComponents[i]!=o.mComponents[i])
                return this->mComponents[i]>o.mComponents[i];
        }
        return false;
    }
    bool operator<=(Vec const&o)const
    {
        if((*this)==o)
            return true;
        
        for(unsigned i=0;i<this->mComponentCount;++i)
        {
            if(this->mComponents[i]!=o.mComponents[i])
                return this->mComponents[i]<o.mComponents[i];
        }
        return false;
    }
    bool operator>=(Vec const&o)const
    {
        if((*this)==o)
            return true;
        
        for(unsigned i=0;i<this->mComponentCount;++i)
        {
            if(this->mComponents[i]!=o.mComponents[i])
                return this->mComponents[i]>o.mComponents[i];
        }
        return false;
    }
    
    float InteriorAngleWith(Vec const&o)const
    {
        assert((this->mComponentCount==2)&&(o.mComponentCount==2));
        
        if((this->Length() == 0) || (o.Length() == 0))
            return 0.0f;
        
        Vec<VecBase2<float> > thisNormalized = this->Normalized();
        Vec<VecBase2<float> > oNormalized = o.Normalized();
        
        float thisAngle = ::atan2(oNormalized.y, oNormalized.x);
        float oAngle = ::atan2(thisNormalized.y, thisNormalized.x);
        
        float retAngle = ::abs(thisAngle - oAngle);
        
        if(retAngle > PI)
            return (2*PI) - retAngle;
        
        return retAngle;
    }
    
    Vec Cross(Vec const&o)
    {
        assert((this->mComponentCount==3)&&(o.mComponentCount==3));
        return Vec(	this->y*o.z - this->z*o.y,
                   this->z*o.x - this->x*o.z,
                   this->x*o.y - this->y*o.x);
    }
    
    T Dot(Vec const&o)const
    {
        T accumulator=0;
        for(unsigned i=0;i<this->mComponentCount;++i)
            accumulator+=this->mComponents[i]*o.mComponents[i];
        return accumulator;
    }
    
    T &operator[](unsigned index)
    {
        assert(index<this->mComponentCount);
        return this->mComponents[index];
    }
    
    T const&operator[](unsigned index)const
    {
        assert(index<this->mComponentCount);
        return this->mComponents[index];
    }
    
    Vec operator-()const
    {
        Vec ret;
        for(unsigned i=0;i<this->mComponentCount;++i)
            ret[i]=-(*this)[i];
        return ret;
    }
    
    Vec &operator+=(Vec const&o)
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
            this->mComponents[i]+=o.mComponents[i];
        
        return *this;
    }
    
    Vec &operator*=(Vec const&o)
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
            this->mComponents[i]*=o.mComponents[i];
        
        return *this;
    }
    
    Vec &operator/=(Vec const&o)
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
            this->mComponents[i]/=o.mComponents[i];
        
        return *this;
    }
    
    Vec &operator-=(Vec const&o)
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
            this->mComponents[i]-=o.mComponents[i];
        
        return *this;
    }
    
    Vec &operator/=(T scalar)
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
            this->mComponents[i]/=scalar;
        
        return *this;
    }
    
    Vec &operator*=(T scalar)
    {
        for(unsigned i=0;i<this->mComponentCount;++i)
            this->mComponents[i]*=scalar;
        
        return *this;
    }
    
    Vec operator+(Vec const&o)const
    {
        Vec ret;
        for(unsigned i=0;i<this->mComponentCount;++i)
            ret[i]=(*this)[i]+o[i];
        return ret;
    }
    
    Vec operator-(Vec const&o)const
    {
        Vec ret;
        for(unsigned i=0;i<this->mComponentCount;++i)
            ret[i]=(*this)[i]-o[i];
        return ret;
    }
    
    /**
     * Yeah, this is a useful operation.
     */
    Vec operator*(Vec const&o)const
    {
        Vec ret;
        for(unsigned i=0;i<this->mComponentCount;++i)
            ret[i]=(*this)[i]*o[i];
        return ret;
    }
    Vec operator/(Vec const&o)const
    {
        Vec ret;
        for(unsigned i=0;i<this->mComponentCount;++i)
            ret[i]=(*this)[i]/o[i];
        return ret;
    }
				
    Vec operator*(typename Base::Type scalar)const
    {
        Vec ret;
        for(unsigned i=0;i<this->mComponentCount;++i)
            ret[i]=(*this)[i]*scalar;
        return ret;
    }
    
    Vec operator/(typename Base::Type scalar)const
    {
        Vec ret;
        for(unsigned i=0;i<this->mComponentCount;++i)
            ret[i]=(*this)[i]/scalar;
        return ret;
    }
    
    T Length()const
    {
        T accum=0;
        for(unsigned i=0;i<this->mComponentCount;++i)
            accum+=(*this)[i]*(*this)[i];
        return ::sqrt(accum);
    }
    
    T SquaredLength()const
    {
        T accum=0;
        for(unsigned i=0;i<this->mComponentCount;++i)
            accum+=(*this)[i]*(*this)[i];
        return accum;
    }
    
    Vec Normalized()const
    {
        return (*this)/Length();
    }
    
    T LargestComponent()const
    {
        assert(this->mComponentCount);
        T ret=this->mComponents[0];
        for(unsigned i=1;i<this->mComponentCount;++i)
        {
            if(this->mComponents[i]>ret)
                ret=this->mComponents[i];
        }
        return ret;
    }
    
    T ComponentWithGreatestMagnitude()const
    {
        assert(this->mComponentCount);
        T ret=this->mComponents[0];
        for(unsigned i=1;i<this->mComponentCount;++i)
        {
            if(abs(this->mComponents[i])>abs(ret))
                ret=this->mComponents[i];
        }
        return ret;
    }
};



template<typename T>
struct VecBase3
{
    typedef T Type;
    static const unsigned mComponentCount=3;
    
    VecBase3()
    {
        for(unsigned i=0;i<mComponentCount;++i)
            mComponents[i]=0;
    }
    
    T				&GetComponent(unsigned index)
    {
        return mComponents[index];
    }
    
    T		   const&GetComponentConst(unsigned index)const
    {
        return mComponents[index];
    }
    
    union
    {
        struct
        {
            T x,y,z;
        };
        struct
        {
            T width,height,depth;
        };
        struct
        {
            T row,column,plane;
        };
        struct
        {
            T red,green,blue;
        };
        T mComponents[3];
    };
};

template<typename T>
struct VecBase4
{
    typedef T Type;
    static const unsigned mComponentCount=4;
    
    VecBase4()
    {
        for(unsigned i=0;i<mComponentCount;++i)
            mComponents[i]=0;
    }
    
    T				&GetComponent(unsigned index)
    {
        return mComponents[index];
    }
    
    T		   const&GetComponentConst(unsigned index)const
    {
        return mComponents[index];
    }
    
    inline Vec<VecBase3<T> > const&xyz()const
    {
        return *reinterpret_cast<Vec<VecBase3<T> > const*>(&mComponents[0]);
    }
    
    union
    {
        struct
        {
            T x,y,z,w;
        };
        struct
        {
            T width,height,depth,wDepth;
        };
        struct
        {
            T row,column,plane,wPlane;
        };
        struct
        {
            T red,green,blue,alpha;
        };
        
        T mComponents[4];
    };
};


template<typename T,unsigned N>
struct VecBaseN
{
    typedef T Type;
    static const unsigned mComponentCount=N;
    
    // Danger: Components are not set to 0
    // TODO: This may be fixable with template specialization (and without "memsetting")
    VecBaseN()
    {
    }
    
    T				&GetComponent(unsigned index)
    {
        return mComponents[index];
    }
    
    T		   const&GetComponentConst(unsigned index)const
    {
        return mComponents[index];
    }
    
    Vec<VecBase3<float> >	GetXYZ()
    {
        assert(mComponents>=3);
        return Vec<VecBase3<float> >(mComponents[0],mComponents[1],mComponents[2]);
    }
    
    T mComponents[N];
};


typedef Vec<VecBaseN<int,1 > >Vec1i;
typedef Vec<VecBase2<int> >	Vec2i;
typedef Vec<VecBase3<int> >	Vec3i;

typedef Vec<VecBaseN<float,1 > >Vec1f;
typedef Vec<VecBase2<float> >	Vec2f;
typedef Vec<VecBase3<float> >	Vec3f;
typedef Vec<VecBase4<float  > >	Vec4f;

typedef Vec<VecBase2<double> >	Vec2d;
typedef Vec<VecBase3<double> >	Vec3d;
typedef Vec<VecBase4<double  > >Vec4d;


template<typename Base>
struct Extrema
{
    Extrema()
    {
    }
    Extrema(Vec<Base> const&min,Vec<Base> const&max)
    : mMin(min), mMax(max)
    {
    }
    Extrema(Extrema const&o)
    : mMin(o.mMin), mMax(o.mMax)
    {
    }
    
    bool operator==(Extrema const&o)
    {
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            if(mMin!=o.mMin)
                return false;
            if(mMax!=o.mMax)
                return false;
        }
        return true;
    }
    
    bool operator!=(Extrema const&o)
    {
        return !((*this)==o);
    }
    
    Extrema operator+(Extrema const&o)const
    {
        return Extrema(mMin + o.mMin, mMax + o.mMax);
    }
    
    Extrema operator*(typename Base::Type scalar)const
    {
        return Extrema(mMin * scalar, mMax * scalar);
    }
    
    /**
     * This function expands the Extrema as needed to enclose a point
     */
    void DoEnclose(Vec<Base> const&pt)
    {
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            if(pt.GetComponentConst(i)<mMin.GetComponentConst(i))
                mMin.GetComponent(i)=pt.GetComponentConst(i);
            
            if(pt.GetComponentConst(i)>mMax.GetComponentConst(i))
                mMax.GetComponent(i)=pt.GetComponentConst(i);
        }
    }
    
    /**
     * This function returns true if the given point is within the region,
     *		false otherwise.
     */
    bool DoesEnclose(Vec<Base> const&pt)const
    {
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            typename Vec<Base>::T minComponent = mMin.GetComponentConst(i);
            typename Vec<Base>::T maxComponent = mMax.GetComponentConst(i);
            typename Vec<Base>::T ptComponent = pt.GetComponentConst(i);
            
            if(ptComponent<minComponent)
                return false;
            
            if(ptComponent>maxComponent)
                return false;
        }
        return true;
    }
    
    /**
     * This function finds the "union" of two sets of extrema, that is the Extrema that enclose both
     */
    Extrema Union(Extrema const&o)const
    {
        Extrema ret;
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            ret.mMin.GetComponent(i)=std::min(mMin.GetComponentConst(i), o.mMin.GetComponentConst(i));
            ret.mMax.GetComponent(i)=std::max(mMax.GetComponentConst(i), o.mMax.GetComponentConst(i));
        }
        return ret;
    }
    
    /**
     * Gets the overlap (region covered by both extrema)
     */
    Extrema GetOverlap(Extrema const&o)const
    {
        Extrema ret;
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            ret.mMin.GetComponent(i)=(mMin.GetComponentConst(i)>o.mMin.GetComponentConst(i))?mMin.GetComponentConst(i):o.mMin.GetComponentConst(i);
            ret.mMax.GetComponent(i)=(mMax.GetComponentConst(i)<o.mMax.GetComponentConst(i))?mMax.GetComponentConst(i):o.mMax.GetComponentConst(i);
            
            if(ret.mMin[i] > ret.mMax[i])
                ret.mMin[i] = ret.mMax[i];
        }
        return ret;
    }
    
    bool GetOverlap(Extrema const&o,Extrema &out)const
    {
        Extrema ret;
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            if( (mMin.GetComponentConst(i) >= o.mMax.GetComponentConst(i)) ||
               (mMax.GetComponentConst(i) <= o.mMin.GetComponentConst(i)))
                return false;
            
            ret.mMin.GetComponent(i)=(mMin.GetComponentConst(i)>o.mMin.GetComponentConst(i))?mMin.GetComponentConst(i):o.mMin.GetComponentConst(i);
            ret.mMax.GetComponent(i)=(mMax.GetComponentConst(i)<o.mMax.GetComponentConst(i))?mMax.GetComponentConst(i):o.mMax.GetComponentConst(i);
        }
        
        out=ret;
        return true;
    }
    
    /**
     * Gets the size (width/height) of the region
     */
    Vec<Base> GetSize()const
    {
        return mMax-mMin;
    }
    
    typename Base::Type GetVolume()const
    {
        Vec<Base> size=GetSize();
        
        typename Base::Type ret = 1;
        
        for(unsigned i=0;i<Base::mComponentCount;++i)
            ret *= size[i];
        
        return ret;
    }
    
    inline bool DoesOverlap(Extrema const&o)const
    {
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            if(mMax.mComponents[i] < o.mMin.mComponents[i])
            {
                return false;
            }
            if(mMin.mComponents[i] > o.mMax.mComponents[i])
            {
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * This function calcuates the corners of the Axially Aligned Bounding Box that encompasses
     *  the extrema.
     *
     * In 1d, the array must be 2 elements, In 2d, 4 elements, in 3d, 8 elements, and so on.
     *
     * NOTE: This function "only" works up to 31 dimensions.
     */
    void CalculateAABB(Vec<Base> *dst)const
    {
        unsigned numDimensions=Base::mComponentCount;
        
        // This variable is a bitfield with one bit per dimension. If that bit is 0,
        //  then the minimum on that dimension is used. If it is 1, then the maximum is used.
        uint32 dimensionMask=0;
        
        // This is a mask with the bit after the last dimension bit set. When we reach this value, we stop the loop.
        uint32 dimensionMaskTooBig=1<<(numDimensions+1);
        
        for(;dimensionMask<dimensionMaskTooBig;++dimensionMask)
        {
            Vec<Base> thisPoint;
            
            for(unsigned currentDimension=0;currentDimension<numDimensions;++currentDimension)
            {
                // This variable is AND'ed with the dimensionMask for each dimension
                uint32 thisDimensionMask=1<<currentDimension;
                
                thisPoint[currentDimension]=
                (dimensionMask&thisDimensionMask)?
                mMax[currentDimension]:mMin[currentDimension];
            }
            
            dst[dimensionMask]=thisPoint;
        }
    }
    
    /**
     * Returns true if the Extrema is valid (minimums <= maximums)
     */
    bool IsValid()
    {
        for(unsigned i=0;i<Base::mComponentCount;++i)
        {
            if(mMin[i] > mMax[i])
                return false;
        }
        
        return true;
    }
    
    Vec<Base>	mMin,mMax;
};

typedef Extrema<VecBaseN<float, 1> >	Extrema1f;
typedef Extrema<VecBaseN<int, 1> >		Extrema1i;

typedef Extrema<VecBase2<int> >			Extrema2i;
typedef Extrema<VecBase3<int> >			Extrema3i;

typedef Extrema<VecBase2<float> >		Extrema2f;
typedef Extrema<VecBase2<double> >		Extrema2d;

typedef Extrema<VecBase3<float> >		Extrema3f;
typedef Extrema<VecBase3<double> >		Extrema3d;

typedef Extrema<VecBase4<float> >		Extrema4f;
typedef Extrema<VecBase4<double> >		Extrema4d;

#endif
