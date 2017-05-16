#ifndef NERDMATRIX_EXCEPTION_H
#define NERDMATRIX_EXCEPTION_H

#include <exception>
#include <string>

namespace nerd {
namespace except {

class NerdMatrixException : public std::exception {
 public:
  NerdMatrixException() : msg_("NerdMatrixException") {}
  explicit NerdMatrixException(const char* msg)
    : msg_(std::string("NerdMatrixException: ") + msg) {}
  explicit NerdMatrixException(const std::string& msg)
    : msg_(std::string("NerdMatrixException: ") + msg) {}
  virtual ~NerdMatrixException() {}
  virtual const char* what() const noexcept {
    return msg_.c_str();
  }

 protected:
  std::string msg_;
};

} /* end of except namespace */
} /* end of nerd namespace */

#endif /* end of include guard: NERDMATRIX_EXCEPTION_H */
