
#include <iostream>
#include <fstream>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include "../src/mat.h"
#include "../src/mat_io.h"

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/shape.hpp>

typedef core::mat<double> mat;
typedef double (*func)(double);
typedef mat (*mfunc)(mat,mat);
typedef mat (*mfuncd)(mat,mat,mat,mat);
using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::string;
using std::vector;

mat addcol(const mat m, const double val) {
  mat newmat(m.nrows, m.ncols+1);
  for (uint32_t i = 0 ; i < m.nrows ; ++i) {
    newmat.data[i * newmat.ncols] = val;
    for (uint32_t j = 0 ; j < m.ncols ; ++j) {
      newmat.data[i * newmat.ncols + j + 1] = m.data[i * m.ncols + j];
    }
  }
  return newmat;
}
mat rmcol(const mat m) {
  mat newmat(m.nrows, m.ncols-1);
  for (uint32_t i = 0 ; i < newmat.nrows ; ++i) {
    for (uint32_t j = 0 ; j < newmat.ncols ; ++j) {
      newmat.data[i * newmat.ncols + j] = m.data[i * m.ncols + j + 1];
    }
  }
  return newmat;
}
class Layer {
 public:
  Layer();
  Layer(const Layer &l);
  Layer(const int pnnodes, const int nnodes, const double lrate,
        const double lambda, double (*act)(double), double (*actd)(double));
  virtual void operator=(const Layer &l);
  virtual mat forwardprop(const mat pa);
  virtual mat backprop(const mat delta);
  virtual void update();
  int getpnnodes() const;
  double getlrate() const;
  double getlambda() const;
  mat getz() const;
  mat geta() const;
  mat getw() const;
  mat getgrad() const;
  mat getdelta() const;
  func getact() const;
  func getactd() const;
  void setw(mat w);

 protected:
  int pnnodes;
  double lrate;
  double lambda;
  double (*act)(double);
  double (*actd)(double);
  mat pa, z, a;
  mat W, grad, delta;
};
class InputLayer : public Layer {
 public:
  InputLayer();
  InputLayer(const InputLayer &input);
  InputLayer(const int innodes);
  virtual void operator=(const InputLayer &input);
  virtual mat forwardprop(const mat input);
};
class OutputLayer : public Layer {
 public:
  OutputLayer();
  OutputLayer(const OutputLayer &output);
  OutputLayer(const int nnodes, const int outputnodes, const double lrate,
              double (*act)(double),
              double (*actd)(double),
              mat (*cost)(mat,mat),
              mat (*costd)(mat,mat,mat,mat));
  virtual void operator=(const OutputLayer &output);
  virtual mat backprop(const mat label);
  mat argmax() const;
  double getcostval() const;
  mfunc getcost() const;
  mfuncd getcostd() const;

 private:
  mat (*cost)(mat,mat);
  mat (*costd)(mat,mat,mat,mat);
  mat y;
};
class NeuralNet {
 public:
  NeuralNet(const InputLayer &input, const OutputLayer &output,
            std::vector<Layer> &layers);
  void feeddata(const mat x, const mat y, const bool check);
  mat predict(const mat sample);
  void gradcheck();
  double computecost();
  mat getresult() const;

 private:
  double computecost(const mat perturb, const uint32_t idx);
  mat computengrad(const int nrows, const int ncols, const int idx);

  const double eps;
  mat x, y;
  mat result;
  mat (*cost)(mat,mat);
  mat (*costd)(mat,mat,mat,mat);
  InputLayer input;
  std::vector<Layer> hidden;
  OutputLayer output;
};
Layer::Layer() : pnnodes(0), lrate(0), lambda(0) {
  act = [] (double x) {return x;};
  actd = [] (double x) {return (x=1);};
}
Layer::Layer(const Layer &l) :
    pnnodes(l.getpnnodes()), lrate(l.getlrate()), lambda(l.getlambda()) {
  z = l.getz();
  a = l.geta();
  W = l.getw();
  grad = l.getgrad();
  act = l.getact();
  actd = l.getactd();
}
Layer::Layer(const int pnnodes, const int nnodes, const double lrate,
             const double lambda, double (*act)(double), double (*actd)(double)) :
             pnnodes(pnnodes+1), lrate(lrate), lambda(lambda), act(act), actd(actd) {
  W = mat(this->pnnodes, nnodes);
  grad = mat(this->pnnodes, nnodes);
  for (uint32_t i = 0 ; i < W.nrows * W.ncols ; ++i) {
    W.data[i] = rand() % 10000 / 10000.0;
  }
}
void Layer::operator= (const Layer &l) {
  pnnodes = l.getpnnodes();
  lrate = l.getlrate();
  lambda = l.getlambda();
  act = l.getact();
  actd = l.getactd();
  z = l.getz();
  a = l.geta();
  W = l.getw();
  grad = l.getgrad();
  delta = l.getdelta();
}
mat Layer::forwardprop(const mat pa) {
  this->pa = addcol(pa, 1);
  z = this->pa * W;
  a = z.f(act);
  return a;
}
mat Layer::backprop(const mat d) {
  // compute this delta and grad
  mat actdz = z.f(actd);
  mat delta = d;
  delta = delta % addcol(actdz, 1);

  delta = rmcol(delta);
  grad = pa.t() * delta;
  // regularization
  grad = grad + W * lambda;

  // compute new delta to throw to next layer
  mat newdelta = delta * W.t();
  return newdelta;
}
void Layer::update() {
  W = W - grad * lrate;
}
int Layer::getpnnodes() const {
  return pnnodes;
}
double Layer::getlrate() const {
  return lrate;
}
double Layer::getlambda() const {
  return lambda;
}
mat Layer::getz() const {
  return z;
}
mat Layer::geta() const {
  return a;
}
mat Layer::getw() const {
  return W;
}
mat Layer::getgrad() const {
  return grad;
}
mat Layer::getdelta() const {
  return delta;
}
func Layer::getact() const {
  return act;
}
func Layer::getactd() const {
  return actd;
}
void Layer::setw(const mat w) {
  this->W = w;
}
InputLayer::InputLayer() {}
InputLayer::InputLayer(const InputLayer &input) :
    Layer(input.getpnnodes(), input.getw().ncols, input.getlrate(),
          input.getlambda(), input.getact(), input.getactd()) {
}
InputLayer::InputLayer(const int innodes) {
  pnnodes = innodes;
  lrate = 1;
  this->act = [](double x) {return x;};
  this->actd = [](double x) {return (x=1);};
  // not using W and grad for input layer
}
void InputLayer::operator= (const InputLayer &input) {
  pnnodes = input.getpnnodes();
  lrate = input.getlrate();
  lambda = input.getlambda();
  act = input.getact();
  actd = input.getactd();
  z = input.getz();
  a = input.geta();
  W = input.getw();
  grad = input.getgrad();
  delta = input.getdelta();
}
mat InputLayer::forwardprop(const mat input) {
  z = input;
  a = input;
  return a;
}
OutputLayer::OutputLayer() {}
OutputLayer::OutputLayer(const OutputLayer &output) {
  cout << "output & constructor" << endl;
  pnnodes = output.getpnnodes();
  lrate = output.getlrate();
  act = output.getact();
  actd = output.getactd();
  z = output.getz();
  a = output.geta();
  W = output.getw();
  grad = output.getgrad();
  delta = output.getdelta();
  cost = output.getcost();
  costd = output.getcostd();
}
OutputLayer::OutputLayer(const int pnnodes, const int outputnodes,
                         const double lrate,
                         double (*act)(double),
                         double (*actd)(double),
                         mat (*cost)(mat,mat),
                         mat (*costd)(mat,mat,mat,mat)) {
  this->pnnodes = pnnodes + 1;
  this->lrate = lrate;
  this->act = act;
  this->actd = actd;
  this->cost = cost;
  this->costd = costd;

  W = mat(this->pnnodes, outputnodes);
  grad = mat(this->pnnodes, outputnodes);
  for (uint32_t i = 0 ; i < W.nrows * W.ncols; ++i) {
    W.data[i] = rand() % 10000 / 10000.0;
  }
}
void OutputLayer::operator= (const OutputLayer &output) {
  pnnodes = output.getpnnodes();
  lrate = output.getlrate();
  lambda = output.getlambda();
  act = output.getact();
  actd = output.getactd();
  z = output.getz();
  a = output.geta();
  W = output.getw();
  grad = output.getgrad();
  delta = output.getdelta();
  cost = output.getcost();
  costd = output.getcostd();
}
mat OutputLayer::backprop(const mat label) {
  y = label;
  delta = costd(y, a, z, pa);
  grad = pa.t() * delta;
  // regularization
  grad = grad + W * lambda;
  // compute next layer delta
  mat newdelta = delta * W.t();
  return newdelta;
}
mat OutputLayer::argmax() const {
  mat result(a.nrows, 2);
  for (uint32_t i = 0 ; i < a.nrows ; ++i) {
    double maxval = a.data[i * a.ncols];
    int maxidx = 0;
    for (uint32_t j = 1 ; j < a.ncols ; ++j) {
      if (a(i, j) > maxval) {
        maxval = a(i, j);
        maxidx = j;
      }
    }
    result.data[i * result.ncols] = maxidx;
    result.data[i * result.ncols + 1] = maxval;
  }
  return result;
}
double OutputLayer::getcostval() const {
  mat J = cost(y, a);
  double val = 0;
  for (uint32_t i = 0 ; i < J.nrows ; ++i) {
    for (uint32_t j = 0 ; j < J.ncols ; ++j) {
      val += J(i, j);
    }
  }
  return val;
}
mfunc OutputLayer::getcost() const {
  return cost;
}
mfuncd OutputLayer::getcostd() const {
  return costd;
}
NeuralNet::NeuralNet(const InputLayer &input, const OutputLayer &output,
                     vector<Layer> &hidden) :
                     eps(1e-5), cost(output.getcost()), costd(output.getcostd()),
                     input(input), hidden(hidden), output(output) {
  for (uint32_t i = 0 ; i < hidden.size() ; ++i)
    cout << hidden[i].getw() << endl;;
}
void NeuralNet::feeddata(const mat x, const mat y, const bool check) {
  this->x = x;
  this->y = y;

  // forward propagation
  mat current = input.forwardprop(x);
  for (uint32_t i = 0 ; i < hidden.size() ; ++i) {
    mat n = hidden[i].forwardprop(current);
    current = n;
  }
  result = output.forwardprop(current);

  // backpropagation
  mat currentdelta = output.backprop(y);
  for (int i = hidden.size()-1 ; i >= 0 ; --i) {
    mat p = hidden[i].backprop(currentdelta);
    currentdelta = p;
  }

  if (check) gradcheck();

  // update parameters
  output.update();
  for (uint32_t i = 0 ; i < hidden.size() ; ++i) {
    hidden[i].update();
  }
}
mat NeuralNet::predict(const mat sample) {
  // prediction is simply forwardprop
  mat current = input.forwardprop(sample);
  for (uint32_t i = 0 ; i < hidden.size() ; ++i) {
    mat n = hidden[i].forwardprop(current);
    current = n;
  }
  result = output.forwardprop(current);
  mat argmax = output.argmax();
  return argmax;
}
void NeuralNet::gradcheck() {
  // back prop result from output to input
#ifdef DEBUG
  const double diff = 1e-6;
  mat w = output.getw();
  mat ngrad = computengrad(w.nrows, w.ncols, hidden.size());
  mat grad = output.getgrad();
  cout << "gradient checking ......";
  for (uint32_t i = 0 ; i < grad.nrows ; ++i) {
    for (uint32_t j = 0 ; j < grad.ncols ; ++j) {
      if (ngrad(i, j) > grad(i, j)) {
        if (ngrad(i, j) - grad(i, j) >= diff) {
          cout << " failed" << endl;
          cout << grad;
          cout << ngrad;
          return;
        }
      } else {
        if (grad(i, j) - ngrad(i, j) >= diff) {
          cout << " failed" << endl;
          cout << grad;
          cout << ngrad;
          return;
        }
      }
    }
  }
  for (int i = hidden.size()-1 ; i >= 0 ; --i) {
    w = hidden[i].getw();
    ngrad = computengrad(w.nrows, w.ncols, i);
    grad = hidden[i].getgrad();

    for (uint32_t i = 0 ; i < grad.nrows ; ++i) {
      for (uint32_t j = 0 ; j < grad.ncols ; ++j) {
        if (ngrad(i, j) > grad(i, j)) {
          if (ngrad(i, j) - grad(i, j) >= diff) {
            cout << " failed" << endl;
            cout << grad;
            cout << ngrad;
            return;
          }
        } else {
          if (grad(i, j) - ngrad(i, j) >= diff) {
            cout << " failed" << endl;
            cout << grad;
            cout << ngrad;
            return;
          }
        }
      }
    }
  }
  cout << " passed." << endl;
#endif
}
double NeuralNet::computecost() {
  mat J = cost(y, result);
  double val = 0;
  for (uint32_t i = 0 ; i < J.nrows * J.ncols ; ++i) {
    val += J.data[i];
  }
  return val;
}
mat NeuralNet::getresult() const {
  return result;
}
double NeuralNet::computecost(const mat perturb, const uint32_t idx) {
  if (idx > hidden.size()) return DBL_MAX;

  // forward propagation
  mat current = input.forwardprop(x);
  for (uint32_t i = 0 ; i < hidden.size() ; ++i) {
    if (i == idx) {
      mat tempw = hidden[i].getw();

      hidden[i].setw(tempw + perturb);
      mat n = hidden[i].forwardprop(current);
      current = n;
      hidden[i].setw(tempw);
    } else {
      mat n = hidden[i].forwardprop(current);
      current = n;
    }
  }
  mat out;
  if (idx == hidden.size()) {
    mat tempw = output.getw();
    output.setw(tempw + perturb);
    out = output.forwardprop(current);
    output.setw(tempw);
  } else {
    out = output.forwardprop(current);
  }

  mat J = cost(y, out);
  double val = 0;
  for (uint32_t i = 0 ; i < J.nrows ; ++i) {
    for (uint32_t j = 0 ; j < J.ncols ; ++j) {
      val += J(i, j);
    }
  }

  // regularization
  for (uint32_t i = 0 ; i < hidden.size() ; ++i) {
    mat tempw = hidden[i].getw();
    if (i == idx) {
      tempw = tempw + perturb;
    }

    const mat regterm =  (tempw % tempw) * (hidden[i].getlambda() / 2.0);
    for (uint32_t j = 0 ; j < regterm.nrows ; ++j) {
      for (uint32_t k = 0 ; k < regterm.ncols ; ++k) {
        val += regterm(j, k);
      }
    }
  }
  return val;
}
mat NeuralNet::computengrad(const int nrows, const int ncols, const int idx) {
  mat wgrad = mat::zero(nrows, ncols);
  mat perturb = mat::zero(nrows, ncols);

  for (int i = 0 ; i < nrows ; ++i) {
    for (int j = 0 ; j < ncols ; ++j) {
      perturb.data[i * perturb.ncols + j] = eps;
      const double loss1 = computecost(perturb, idx);
      const double loss2 = computecost(-perturb, idx);
      wgrad.data[i * wgrad.ncols + j] = (loss1 - loss2) / (2 * eps);

      perturb.data[i * perturb.ncols + j] = 0;
    }
  }

  return wgrad;
}

double rectifier(double z) {
  return z >= 0 ? z : 0;
}
double rectifiergrad(double z) {
  return z >= 0 ? 1 : 0;
}
double sigmoid(double z) {
  return 1.0 / (1.0 + exp(-z));
}
double sigmoidgrad(double z) {
  const double e = exp(-z);
  const double b = 1 + e;
  return e / (b * b);
}
mat cost(mat y, mat h) {
  mat one = mat::ones(y.nrows, y.ncols);
  mat J = -(y % h.f(log) + (one-y) % (one-h).f(log));
  return J;
}
mat costd(mat y, mat a, mat,mat) {
  mat grad = (a - y);
  return grad;
}
void load(mat &x, mat &y) {
  std::ifstream xinput("/home/joseph/C/project/nn/data/samplex.data", std::ios::in);
  std::ifstream yinput("/home/joseph/C/project/nn//data/sampley.data", std::ios::in);
  x = mat(100, 2);
  y = mat::zero(100, 2);
  if (xinput.is_open() && yinput.is_open()) {
    cout << "reading data..." << endl;
    for (uint32_t i = 0 ; i < 100 ; i ++) {
      double x0 = 0, x1 = 0, x2 = 0;
      int yi = 0;
      xinput >> x0;
      xinput >> x1;
      xinput >> x2;
      yinput >> yi;
      x.data[i * x.ncols] = x1;
      x.data[i * x.ncols + 1] = x2;
      y.data[i * y.ncols + yi] = 1;
    }
  }
}
void loadsample(mat &sample, const int w, const int h) {
  sample = mat(w * h, 2);
  for (int i = 0 ; i < h * w ; ++i) {
    sample.data[i * sample.ncols] = static_cast<double>(1.0 * (i % w) / w);
    sample.data[i * sample.ncols + 1] = static_cast<double>(1.0 * (i / w) / h);
  }
}

int main(void) {
  const double lrate = 1e-3;
  const double lambda = 1e-2;
  const int w = 800;
  const int h = 600;

  srand(time(NULL));

  InputLayer input(2);
  vector<Layer> hidden = {
    //Layer(2, 2, lrate, rectifier, rectifiergrad),
    //Layer(2, 6, lrate, atan, [](double x) {return 1.0/(1.0+x*x);}),
    Layer(2, 6, lrate, lambda, atan, [](double x) {return 1.0/(1.0+x*x);}),
    Layer(6, 6, lrate, lambda, rectifier, rectifiergrad),
    Layer(6, 6, lrate, lambda, sigmoid, sigmoidgrad),
    Layer(6, 6, lrate, lambda, atan, [](double x) {return 1.0/(1.0+x*x);}),
    //Layer(2, 6, lrate, sigmoid, sigmoidgrad),
    //Layer(6, 6, lrate, sigmoid, sigmoidgrad),
    //Layer(6, 6, lrate, sigmoid, sigmoidgrad),
    //Layer(6, 6, lrate, sigmoid, sigmoidgrad),
  };
  OutputLayer output(6, 2, lrate, sigmoid, sigmoidgrad, cost, costd);
  NeuralNet nnet(input, output, hidden);

  mat x, y, sample;
  load(x, y); loadsample(sample, w, h);
  nnet.feeddata(x, y, true);
  for (int i = 0 ; i < 120000 ; ++i) {
    nnet.feeddata(x, y, false);
    //nnet.feeddata(x, y, true);
    cout << "\riteration: " << i << " | cost: " << nnet.computecost();
  }
  cout << endl;
  mat result = nnet.predict(sample);

  cv::Mat canvas = cv::Mat::zeros(h, w, CV_8UC1);
  for (uint32_t i = 0 ; i < result.nrows ; ++i) {
    canvas.at<uchar>(i/w, i%w) = result(i, 0) == 1 ? 128 : 0;
  }
  for (uint32_t i = 0 ; i < x.nrows ; ++i) {
    if (y(i,0) == 0) {
      cv::circle(canvas, cv::Point(x(i,0)*w, x(i,1)*h), 3, cv::Scalar(255));
    } else {
      cv::circle(canvas, cv::Point(x(i,0)*w, x(i,1)*h), 3, cv::Scalar(100));
    }
  }

  cv::imshow("demo", canvas);
  cv::waitKey(0);

  return 0;
}
