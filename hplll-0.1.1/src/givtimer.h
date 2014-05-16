/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* kernel/system/givtimer.h
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added _start_t member to BaseTimer, so that stop () does not clobber the
 * class' memory of its start time. This allows it to be called repeatedly to
 * get elapsed times.
 * ------------------------------------
 *
 * See COPYING for license information.
 *
 */

/** @file givtimer.h
 * @ingroup system
 * @brief timer
 */

#ifndef __GIVARO_timer_H
#define __GIVARO_timer_H

//#include <iostream>
//#include <givaro/givconfig.h>

namespace Givaro {
// class BaseTimer; class RealTimer; class SysTimer; class UserTimer;

/** \brief base for class RealTimer; class SysTimer; class UserTimer;
  \ingroup util
  */
class BaseTimer {
public:
	enum {
		MSPSEC = 1000000  // microsecond per second
	};

	BaseTimer() : _start_t(0.), _t(0.) {}

	// -- Clear timer :
	inline void clear()
	{
		_t = 0.0;
	}

	// -- total amount of second spent
	inline double time() const
	{
		return _t;
	}

	// -- Return a value to initialize random generator
	static long seed();

	// -- basic methods:
	std::ostream& print( std::ostream& ) const;

	// -- Some arithmetic operators to compute cumulative time :
	BaseTimer& operator = (const BaseTimer & T) ;
	const BaseTimer operator - (const BaseTimer & T)  const;
	const BaseTimer operator - () ;
	const BaseTimer operator +  (const BaseTimer & T)  const;
	BaseTimer& operator += (const BaseTimer & T) { return *this = *this + T; };
	BaseTimer& operator -= (const BaseTimer & T) { return *this = *this - T; };
	public:
	double _start_t;  // time as of start ()
	double _t;        // time
};

//! I/O
inline std::ostream& operator<< (std::ostream& o, const BaseTimer& BT)
{
	return BT.print(o);
}

//! Real timer
class RealTimer : public BaseTimer {
public:
	inline RealTimer( const BaseTimer& BT ): BaseTimer(BT) {};
	inline RealTimer(){};
	void start();
	void stop();
	inline double elapsed()
	{
		stop();
		return _t ; // time()
	}

};


//! User timer
class UserTimer : public BaseTimer {
public:
	inline UserTimer( const BaseTimer& BT ) : BaseTimer(BT) {};
	inline UserTimer() {};
	void start();
	void stop();
	inline double elapsed()
	{
		stop();
		return _t ; // time()
	}

};

//! Sys timer
class SysTimer : public BaseTimer {
public:
	inline SysTimer( const BaseTimer& BT ): BaseTimer(BT) {};
	inline SysTimer() {};
	void start();
	void stop();
	inline double elapsed()
	{
		stop();
		return _t ; // time()
	}

};


//! Timer
class Timer {
public :

	Timer() : _count(0)
	{
		rt.clear();
		ut.clear();
		st.clear();
	}

	/*! Clear timer.
	 * Everything reset to 0. This need not be called before the first
	 * start since the constructor does it.
	 */
	void clear();

	/*! Start timer.
	 * Starts the timer.
	 * If called after  another \c start() or a \c stop(), it sets the timer
	 * to a totally fresh new start.
	 */
	void start();

	/*! Stop timer.
	 * Stops the timer. The time since the previous \c start() is stored.
	 * If called again, \c stop() will store the time since the previous \c
	 * start() again, acting as a \c pause().
	 * @pre \c start() should have been called before...
	 */
	void stop();

	/*! total amount of second spent in user mode.
	 * @return the user time elapsed between the latest \c start() and the
	 * latest \c stop().
	 * @pre \c stop() is called before.
	 */
	double usertime() const { return ut.time(); }

	/*! total amount of second spent in system mode.
	 * @return the system time elapsed between the latest \c start() and
	 * the latest \c stop().
	 * @pre \c stop() is called before.
	 */
	double systime () const { return st.time(); }

	/*! real total amount of second spent.
	 * @return the real total time elapsed between the latest \c start()
	 * and the latest \c stop().
	 * @pre \c stop() is called before.
	 */
	double realtime () const { return rt.time(); }

	/*! User mode time spent since start.
	 * A call to \c stop() is useless.
	 * @return elpased time (in seconds) since \c start() in user mode.
	 */
	double userElapsedTime()  { return ut.elapsed(); }

	/*! System mode time spent since start.
	 * A call to \c stop() is useless.
	 * @return elpased time (in seconds) since \c start() in system mode.
	 */
	double sysElapsedTime ()  { return st.elapsed(); }

	/*! real total amount of second spent since start.
	 * A call to \c stop() is useless.
	 * @return elpased time (in seconds) since \c start().
	 */
	double realElapsedTime () { return rt.elapsed(); }


	// retourne une petite graine
	// long seed() const { return RealTimer::seed(); }

	// Some arithmetic operators to compute cumulative time :
	Timer& operator = (const Timer & T) ;
	const Timer operator - (const Timer & T)  const;
	const Timer operator - () ;
	const Timer operator + (const Timer & T)  const;
	Timer& operator += (const Timer & T) { return *this = *this + T; };
	Timer& operator -= (const Timer & T) { return *this = *this - T; };


	// -- methods :
	std::ostream& print( std::ostream& ) const;

	size_t count() const
	{
		return _count;
	}

private:
	size_t _count; // how many
	RealTimer rt;
	UserTimer ut;
	SysTimer  st;
};

//! I/O
inline std::ostream &operator << (std::ostream &o, const Timer &T)
{
	double ut = T.usertime();
	if (ut < 0.0000000001) ut = 0;
	return o << T.realtime() << "s (" << ut << " cpu) [" << T.count() << "]";
}
}

//#ifdef GIVARO_USES_OMP
#include <omp.h>

namespace Givaro {
//! OMP timer
struct OMPTimer {
	double _c;
	void start() { _c = omp_get_wtime(); }
	void stop() { _c = omp_get_wtime() - _c; }
	void clear() { _c = 0.0; }
	double realtime() { return _c; }
	double usertime() { return _c; }
	double time() const { return _c; }
	friend std::ostream& operator<<(std::ostream& o, const OMPTimer& t) {
		return o << t._c << 's';
	}

	OMPTimer& operator =(const OMPTimer& t) { _c = t._c; return *this; }
	OMPTimer& operator+=(const OMPTimer& t) { _c += t._c; return *this; }
	OMPTimer& operator-=(const OMPTimer& t) { _c -= t._c; return *this; }
	OMPTimer  operator +(const OMPTimer& t) const
	{
		OMPTimer r; r._c = _c + t._c; return r;
	}
	OMPTimer  operator -(const OMPTimer& t) const
	{
		OMPTimer r; r._c = _c - t._c; return r;
	}
	OMPTimer  operator -() { OMPTimer r; r._c = - _c; return r; }
};
} // namespace Givaro
//#endif // OMP


#endif // __GIVARO_timer_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
