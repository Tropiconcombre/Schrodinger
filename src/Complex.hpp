#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <iostream>

class Complex {
  // allows this operator to access the private fields even if it's not part of the class
  friend std::ostream& operator<<(std::ostream& out, const Complex& c);

private:
  const double EPSILON = 1e-4; // The accuracy, to be used when testing equality
  const double r_      = 0;    // real part, 0 by default
  const double i_      = 0;    // imaginary part, 0 by default

public:
  Complex();                   // Constructor with no arguments (creates 0 + 0i)
  Complex(double r, double i); // Constructor

  double real() const; // Returns the real component
  double imag() const; // Returns the imaginary component

  Complex operator+(Complex const& c2) const; // Returns a new complex being the addition of `this` and `c2`
  Complex operator-(Complex const& c2) const;
  Complex operator*(Complex const& c2) const;
  Complex operator/(Complex const& c2) const;
  bool operator==(const Complex& c2) const;
  bool operator!=(const Complex& c2) const;
};

#endif
