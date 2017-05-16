#ifndef NERDMATRIX_ILLEGAL_OPS_H
#define NERDMATRIX_ILLEGAL_OPS_H

#include "NerdMatrix/except/nerdmatrix_exception.h"
#include <exception>
#include <string>

namespace nerd {
namespace except {

class IllegalOpsException : public NerdMatrixException {
 public:
  IllegalOpsException() : msg_("IllegalOpsException") {}
  explicit IllegalOpsException(const char* msg)
      : msg_(std::string("IllegalOpsException: ") + msg) {}
  explicit IllegalOpsException(const std::string& msg)
      : msg_(std::string("IllegalOpsException: ") + msg) {}
  virtual ~IllegalOpsException() {}
  virtual const char* what() const noexcept { return msg_.c_str(); }

 protected:
  std::string msg_;
};

} /* end of except namespace */
} /* end of nerd namespace */

#endif /* end of include guard: NERDMATRIX_ILLEGAL_OPS_H */
