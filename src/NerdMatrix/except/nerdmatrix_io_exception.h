#ifndef NERDMATRIX_IO_EXCEPTION_H
#define NERDMATRIX_IO_EXCEPTION_H

#include "NerdMatrix/except/nerdmatrix_exception.h"

namespace nerd {
namespace except {

class NerdMatrixIOException : public NerdMatrixException {
 public:
  NerdMatrixIOException() : msg_("IllegalOpsException") {}
  explicit NerdMatrixIOException(const char* msg)
      : msg_(std::string("IllegalOpsException: ") + msg) {}
  explicit NerdMatrixIOException(const std::string& msg)
      : msg_(std::string("IllegalOpsException: ") + msg) {}
  virtual ~NerdMatrixIOException() {}
  virtual const char* what() const noexcept { return msg_.c_str(); }

 protected:
  std::string msg_;
};

} /* end of except namespace */
} /* end of nerd namespace */

#endif /* end of include guard: NERDMATRIX_IO_EXCEPTION_H */
