#ifndef __UTL_H_INCLUDED
#define __UTL_H_INCLUDED


/*===========================================================================*/
/* Linear Algebra utilities
/*---------------------------------------------------------------------------*/

#include <Eigen/Dense>
#include <math.h>


namespace la
{

  template<typename T,unsigned int N>
  struct vec_t
  {
    typedef Eigen::Matrix<T,1,N> type;
  };

  typedef vec_t<double,2>::type dvec2_t;
  typedef vec_t<double,3>::type dvec3_t;
  typedef vec_t<double,4>::type dvec4_t;

  typedef vec_t<float,2>::type fvec2_t;
  typedef vec_t<float,3>::type fvec3_t;
  typedef vec_t<float,4>::type fvec4_t;

  typedef vec_t<int,2>::type ivec2_t;
  typedef vec_t<int,3>::type ivec3_t;
  typedef vec_t<int,4>::type ivec4_t;

  typedef vec_t<unsigned int,2>::type uivec2_t;
  typedef vec_t<unsigned int,3>::type uivec3_t;
  typedef vec_t<unsigned int,4>::type uivec4_t;

  template<typename T>
  inline typename vec_t<T,2>::type
  make_vec(const T&a , const T&b )
  {
    typename vec_t<T,2>::type v;
    v[0] = a; v[1] = b; return v;
  }

  template<typename T>
  inline typename vec_t<T,3>::type
  make_vec(const T&a , const T&b ,const T&c)
  {
    typename vec_t<T,3>::type v;
    v[0] = a; v[1] = b; v[2] = c; return v;
  }

  template<typename T>
  inline typename vec_t<T,4>::type
  make_vec(const T&a , const T&b ,const T&c,const T&d)
  {
    typename vec_t<T,4>::type v;
    v[0] = a; v[1] = b; v[2] = c; v[3] = d; return v;
  }

  template<typename T>
  inline typename vec_t<T,4>::type
  make_vec(const typename vec_t<T,3>::type &u, const T&d)
  {
    typename vec_t<T,4>::type v;
    v[0] = u[0]; v[1] = u[1]; v[2] = u[2]; v[3] = d; return v;
  }

  template<typename T>
  inline typename vec_t<T,4>::type
  make_vec(const T&a, const typename vec_t<T,3>::type &u)
  {
    typename vec_t<T,4>::type v;
    v[0] = a; v[1] = u[0]; v[2] = u[1]; v[3] = u[2]; return v;
  }


  template<typename T,unsigned int N>
  inline double tri_area(const typename vec_t<T,N>::type &p,
                    const typename vec_t<T,N>::type &q,
                    const typename vec_t<T,N>::type &r)
  {
    double a = (p-q).norm();
    double b = (q-r).norm();
    double c = (r-p).norm();

    double s = (a+b+c)/2.0;

    double area2  = s*(s-a)*(s-b)*(s-c);

    double area =  sqrt(area2);

    return area;
  }
}

/*===========================================================================*/




/*===========================================================================*/
/* Misc utility functions
/*---------------------------------------------------------------------------*/

namespace utl {

/*---------------------------------------------------------------------------*/

/**
  \brief A generic to string impl that uses the << operator of type T

  \tparam     T  the value type of the vector
  \param[in]  v  Value

  \returns string
**/

template <typename T>
inline std::string to_string (const T& v)
{std::stringstream ss; ss << v ; return ss.str();}

/*---------------------------------------------------------------------------*/

template <class iter_t>
void argsort(iter_t b, iter_t e, std::vector<size_t>& idxs);

/*---------------------------------------------------------------------------*/

template <class iter_t>
void rearrange(iter_t b, iter_t e,const std::vector<size_t>& idxs);

/*---------------------------------------------------------------------------*/

/// \brief strip front and back whitespaces in given string
void trim(std::string &s);

/*---------------------------------------------------------------------------*/

}// namespace utl
/*===========================================================================*/






/*===========================================================================*/
/* Misc utility classes
/*---------------------------------------------------------------------------*/

#include <boost/iterator/iterator_facade.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

namespace utl {

/*---------------------------------------------------------------------------*/


/**
  \brief A class to iterate over non-comment lines of a file

  \note All parts of a line that succeed the comment char ('#' default) are
  assumed to be a comment and are stripped

  \note Sample usage
  for(file_line_iterator lgen("file.txt"),lend; lgen != lend; )
    cout << *lgen++;

**/

class file_line_iterator
    : public boost::iterator_facade<
    file_line_iterator
    , std::string const
    , boost::forward_traversal_tag
    >
{
public:
  file_line_iterator(const char * f,char c_char='#');
  file_line_iterator(){}

private:
  friend class boost::iterator_core_access;
  void increment();
  bool equal(file_line_iterator const& other) const;
  const std::string &dereference() const;

private:
  boost::shared_ptr<std::ifstream> is;
  std::string                      value;
  char                             c_char;
};

/*---------------------------------------------------------------------------*/

/** \brief Simple stopwatch timer (to measure wall clock time NOT cpu time) **/
class timer
{
 public:
  timer()
  {restart();}

  inline void   restart()
  { _start_time = boost::posix_time::microsec_clock::local_time(); }

  inline double elapsed() const
  {
    boost::posix_time::time_duration td =
        boost::posix_time::microsec_clock::local_time() - _start_time;

    return double(td.total_milliseconds())/1000;
  }

 private:
  boost::posix_time::ptime _start_time;
}; // timer

/*---------------------------------------------------------------------------*/

/** \brief A simple multi-threaded logger                                   **/
class logger
{
 public:

  enum severity_level {trace,debug,info,warning,error,fatal};

  inline bool isOpen(severity_level severity)
  {return severity >= info;}

  inline void push(const std::string & log)
  {
    std::string tstr = boost::posix_time::to_simple_string
        (boost::posix_time::second_clock::local_time());

    boost::mutex::scoped_lock  lock(s_mutex);

    std::clog << " [" <<tstr <<"]"
              << " [" <<boost::this_thread::get_id()<<"]"
              << std::endl
              << log
              << std::endl;
  }

  static inline logger& get() {return s_logger;}

private:

  static boost::mutex s_mutex;
  static logger       s_logger;
}; // logger

namespace detail
{
/** \brief a small extension of std::pair where only first is initialized   **/
template<class _T1, class _T2> struct pair:public std::pair<_T1,_T2>
{pair(const _T1 &t1){std::pair<_T1,_T2>::first = t1;}};
}//namespace detail

}// namespace utl

/*===========================================================================*/




/*===========================================================================*/
/* Macro definitions
/*---------------------------------------------------------------------------*/

#define is_in_range(i,b,e) (((b) <= (i)) && ((i) < (e)))

/*---------------------------------------------------------------------------*/

#define LOG_MACRO(lev)  \
  if(utl::logger::get().isOpen(utl::logger::lev))\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push(__utl_lm_v__.second.str())) __utl_lm_v__.second

/*---------------------------------------------------------------------------*/

#define STRMVAR(VAR) #VAR << " = "<< (VAR)

/*---------------------------------------------------------------------------*/

#ifndef ASSERT
#ifndef NDEBUG
#define ASSERT(cond) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERT(cond)
#endif // ifndef NDEBUG
#endif // ifndef ASSERT

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV
#ifndef NDEBUG
#define ASSERTV(cond,var1) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV(cond,var1)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV2
#ifndef NDEBUG
#define ASSERTV2(cond,var1,var2) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV2(cond,var1,var2)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV2

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV3
#ifndef NDEBUG
#define ASSERTV3(cond,var1,var2,var3) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
    << #var3 <<" = "<< (var3) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV3(cond,var1,var2)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV2

/*---------------------------------------------------------------------------*/

#define ENSURE(cond,mesg) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV(cond,mesg,var1) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV2(cond,mesg,var1,var2) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n";\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV3(cond,mesg,var1,var2,var3) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
    << #var3 <<" = "<< (var3) << "\n";\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSURES(cond) \
  if(!(cond)) \
  for(std::stringstream ss ; true ; throw std::runtime_error(ss.str())) \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "

/*===========================================================================*/


#endif
