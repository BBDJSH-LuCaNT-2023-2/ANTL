/**
 * @file timer.c
 * @author Thomas Papanikolaou (TP)
 * @version $Header$
 */

#include <NTL/ZZ.h>
#include "timer.hpp"

namespace ANTL
{

  // Define static member.
  struct rusage timer::buffer;

  void timer::print_hms(ostream & out, long rt) const
  {
    long d, h, m, s, hs;

    if (rt)
      {
	d = rt / 86400000;
	if (d)
	  {
	    out << d << " day ";
	    rt -= d * 86400000;
	  }

	h = rt / 3600000;
	if (h)
	  {
	    out << h << " hour ";
	    rt -= h * 3600000;
	  }

	m = rt / 60000;
	if (m)
	  {
	    out << m << " min ";
	    rt -= m * 60000;
	  }

	s = rt / 1000;
	if (s)
	  {
	    out << s << " sec ";
	    rt -= s * 1000;
	  }

	hs = rt / 10;
	if (hs)
	  {
	    out << hs << " hsec ";
	    rt -= hs * 10;
	  }

	if (rt)
	  out << rt << " msec";
      }
    else
      out << "0 msec";
  }

  timer::timer()
  {
    user_t = 0;
    sys_t = 0;
    t_user_t = 0;
    t_sys_t = 0;
    print_mode = TIME_MODE;
  }

  timer::timer(const timer & t)
  {
    user_t = t.user_t;
    sys_t = t.sys_t;
    t_user_t = t.t_user_t;
    t_sys_t = t.t_sys_t;
    print_mode = t.print_mode;
  }

  timer::~timer()
  {
  }

  int timer::set_print_mode(int m)
  {
    int old_print_mode = print_mode;
    switch (m)
      {
      case 0 :
      case 1 :
	print_mode = m;
	break;
      default:
	print_mode = 1;
	break;
      }
    return old_print_mode;
  }

  int timer::get_print_mode() const
  { return print_mode; }

  void timer::start_timer()
  {
#if defined (__linux__)
    getrusage(RUSAGE_SELF, &timer::buffer);
#else
    getrusage(0, &timer::buffer);
#endif
    user_t = timer::buffer.ru_utime.tv_sec * 100
      + timer::buffer.ru_utime.tv_usec / 10000;
    sys_t  = timer::buffer.ru_stime.tv_sec * 100
      + timer::buffer.ru_stime.tv_usec / 10000;
    t_user_t = 0;
    t_sys_t  = 0;
  }

  void timer::stop_timer()
  {
#if defined (__linux__)
    getrusage(RUSAGE_SELF, &timer::buffer);
#else
    getrusage(0, &timer::buffer);
#endif
    t_user_t += timer::buffer.ru_utime.tv_sec * 100
      + timer::buffer.ru_utime.tv_usec / 10000
      - user_t;
    t_sys_t  += timer::buffer.ru_stime.tv_sec * 100
      + timer::buffer.ru_stime.tv_usec / 10000
      - sys_t;
  }

  void timer::cont_timer()
  {
#if defined (__linux__)
    getrusage(RUSAGE_SELF, &timer::buffer);
#else
    getrusage(0, &timer::buffer);
#endif
    user_t = timer::buffer.ru_utime.tv_sec * 100
      + timer::buffer.ru_utime.tv_usec / 10000;
    sys_t  = timer::buffer.ru_stime.tv_sec * 100
      + timer::buffer.ru_stime.tv_usec / 10000;
  }

  long timer::user_time() const
  { return t_user_t; }

  long timer::sys_time() const
  { return t_sys_t; }

  long timer::real_time() const
  { return t_user_t + t_sys_t; }

  void timer::print(ostream & out) const
  {
    switch (print_mode)
      {
      case 0:
        out << real_time() << " real\t";
        out << user_time() << " user\t";
        out <<  sys_time() << " sys";
	break;
      case 1:
	print_hms(out, real_time() * 10);
	break;
      default:
	print_hms(out, real_time() * 10);
	break;
      }
  }

  timer & timer::operator = (const timer & t)
  {
    user_t = t.user_t;
    sys_t = t.sys_t;
    t_user_t = t.t_user_t;
    t_sys_t = t.t_sys_t;
    print_mode = t.print_mode;
    return *this;
  }

  ostream & operator << (ostream & out, const timer & t)
  { t.print(out); return out; }


  timer operator - (const timer & t1, const timer & t2)
  {
    timer t3;
    t3.t_user_t = t1.t_user_t - t2.t_user_t;
    t3.t_sys_t = t1.t_sys_t - t2.t_sys_t;
    return t3;
  }

  //
  // MyTime
  //
  // Task:
  //      outputs a time in my own special format.
  //

  void
  MyTime(long t)
  {
    MyTime(to<ZZ>(t));
  }

  void
  MyTime(const ZZ & tt)
  {
    ZZ d,h,m,s,ms,t;

    t = tt;
    ms = t % 100;
    t /= 100;
    s = t % 60;
    t /= 60;
    m = t % 60;
    t /= 60;
    h = t % 24;
    d = t / 24;

    if (d > 0) {
      cout << d << " day, ";
      cout << h << " hour, ";
      cout << m << " min, ";
      cout << s << " sec, ";
      cout << ms << " csec";
    }
    else if (h > 0) {
      cout << h << " hour, ";
      cout << m << " min, ";
      cout << s << " sec, ";
      cout << ms << " csec";
    }
    else if (m > 0) {
      cout << m << " min, ";
      cout << s << " sec, ";
      cout << ms << " csec";
    }
    else if (s > 0) {
      cout << s << " sec, ";
      cout << ms << " csec";
    }
    else
      cout << ms << " csec";

    return;
  }

} // ANTL
