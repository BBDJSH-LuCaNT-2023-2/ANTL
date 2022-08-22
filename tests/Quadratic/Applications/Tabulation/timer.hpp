/**
 * @file timer.hpp
 * @author Thomas Papanikolaou (TP)
 * @brief The LiDIA timer.
 * @version $Header$
 * MJJ:  REPLACE THIS WITH TIME IN MSEC STUFF FROM QIR!
 */

#ifndef ANTL_TIMER_H
#define ANTL_TIMER_H

#include <iostream>

#include <sys/time.h>
#include <sys/resource.h>

#include <ANTL/common.hpp>

#define TIME_MODE 0
#define HMS_MODE 1

namespace ANTL
{

#ifdef sun
  extern "C"
  {
    int getrusage (int who, struct rusage *rusage);
  }
#endif

  class timer
  {
  private:

    static struct rusage buffer;

    long user_t, sys_t;
    long t_user_t, t_sys_t;
    int print_mode;

    void print_hms (std::ostream & out = std::cout, long rt = 0) const;

  public:

    timer ();
    timer (const timer &);
    ~timer ();

    int set_print_mode (int x = 1);
    int get_print_mode () const;

    void start_timer ();
    void stop_timer ();
    void cont_timer ();

    long user_time () const;
    long sys_time () const;
    long real_time () const;

    void print (std::ostream & out = std::cout) const;

    timer & operator = (const timer & t);

    friend std::ostream & operator << (std::ostream & out, const timer & t);
    friend timer operator - (const timer & t1, const timer & t2);

    friend void MyTime (long t);
    friend void MyTime (const ZZ & tt);
  };

  void MyTime (long t);
  void MyTime (const ZZ & tt);
  timer operator - (const timer & t1, const timer & t2);

} // ANTL

#endif // guard

