/* 
 * Utility logging file for PIPS-S. Thread-safe. 
 *
 * Caution: There is an overhead associated with logging, therefore avoid excessive logging,
 * especially in critical parts the code. Some overhead occurs even when the logging's filter
 * is set to reduce/disable output.
 *
 * Usage: 
 *
 * Call init_logging(lvl) close to your application entry point. lvl indicates the verbosity level:
 *  0 - debug: high verbosity output, including debug
 *  1 - info: allows informational, application-related messages  and output
 *  2 - summary: less verbose than 'info'. It does NOT include algorithm iteration output and 
 *  application informational logging. Prints a summary of the solve process on exit.
 *  3 - warning: only messages with severity equal to or higher than warning
 *  4 - error: only errors and critical messages
 *  5 - critical: very compact, only critical messages. If execution is successfull, then 
 *  there is no logging.
 *
 * Use one the below  macros for logging. Other BOOST logging macros also will work, however
 * the output may not be formatted properly or some logging message may be lost. 
 * 
 * PIPS_ALG_LOG_SEV(lvl) << msg;
 *  - this macro is intended for ALGORITHM related messages and will simply output 'msg'
 *
 * PIPS_APP_LOG_SEV(lvl) << msg;
 *  - this macro is intended for application messages and will result in a line with the timestamp and
 * severity preceeding the text in 'msg':
 *  [17:19:10.250380] [critical] "critical message contained in msg"
 *
 *
 * Author: Cosmin G. Petra - Argonne National Laboratory, March 2015.
 */

#define BOOST_LOG_DYN_LINK 1

#ifndef PIPSLOGGING
#define PIPSLOGGING

#include <boost/smart_ptr/shared_ptr.hpp>
//#include <boost/core/null_deleter.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/utility/manipulators/add_value.hpp>
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/format.hpp>

namespace logging = boost::log;
namespace expr = boost::log::expressions;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;

namespace boost {

//! A function object that does nothing and can be used as an empty deleter for \c shared_ptr
struct null_deleter
{
    //! Function object result type
    typedef void result_type;
    /*!
     * Does nothing
     */
    template< typename T >
    void operator() (T*) const BOOST_NOEXCEPT {}
};

}

using boost::shared_ptr;


enum severity_level
{
  debug, info, summary,  warning,  error, fatal
};

// The formatting logic for the severity level
template< typename CharT, typename TraitsT >
inline std::basic_ostream< CharT, TraitsT >& operator<< (std::basic_ostream< CharT, TraitsT >& strm, severity_level lvl)
{
  static const char* const str[] = {"debus", "info", "summary", "warning", "error", "critical"};
    if (static_cast< std::size_t >(lvl) < (sizeof(str) / sizeof(*str)))
        strm << str[lvl];
    else
        strm << static_cast< int >(lvl);
    return strm;
}

class PIPSLogging
{
protected:
  //hidden constructors; this class is not to be instantiated
  PIPSLogging() {};
public:
  // severity logger instance
  static src::severity_logger< severity_level > g_sev_log;
  
  //logging functions, however, MACROS defined at the end of the file are faster
  //and easier to use
  static void AppLogSeverity(severity_level lvl, std::string msg) 
  {
    BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl) << msg;
  }
  static void AlgLogSeverity(severity_level lvl, std::string msg)
  {
    BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl) << logging::add_value("ALGORITHM", "on") << msg;
  }

  static void init_logging(int level)
  {
    
    //for now the output will be only on the console
    
    //create a text output sink:
    typedef sinks::synchronous_sink< sinks::text_ostream_backend > text_sink;
    shared_ptr< text_sink > pSink(new text_sink);
    
    // Here synchronous_sink is a sink frontend that performs thread synchronization
    // before passing log records to the backend; this makes backend easier to implement
    
    text_sink::locked_backend_ptr pBackend = pSink->locked_backend();
    
    // Add the stream corresponding to the console
    shared_ptr< std::ostream > pStream(&std::clog, boost::null_deleter());
    pBackend->add_stream(pStream);
    
    // More than one stream to the sink backend. For example to log to a file:
    //shared_ptr< std::ofstream > pStream2(new std::ofstream("sample.log"));
    //assert(pStream2->is_open());
    //pBackend->add_stream(pStream2);
    
    // Add the sink to the logging library
    logging::core::get()->add_sink(pSink);
    //set the formatter for the sink's output
    pSink->set_formatter(expr::stream << 
			 expr::if_(expr::has_attr("ALGORITHM"))
			 [ 
			  expr::stream << "  " << expr::smessage//expr::attr< std::string >("ALGORITHM")
			   ]
			 .else_
			 [
			  expr::stream 
			  << "[" 
			  << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%H:%M:%S.%f") 
			  << "] [" << expr::attr< severity_level >("Severity") << "] "
			  << expr::smessage
			  ]);
    
    pSink->set_filter(expr::attr< severity_level >("Severity").or_default(info) >= level);
    
    //add a time stamp
    attrs::local_clock TimeStamp;
    logging::core::get()->add_global_attribute("TimeStamp", TimeStamp);    
  };
};

#define PIPS_APP_LOG_SEV(lvl) BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl)
#define PIPS_ALG_LOG_SEV(lvl) BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl) << logging::add_value("ALGORITHM", "on")


#endif
