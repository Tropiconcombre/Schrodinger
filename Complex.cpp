#include "Complex.hpp"

Complex::Complex() {
}

Complex::Complex(double r, double i) : r_(r), i_(i) {
}

double Complex::real() const {
   return this->r_;
}
double Complex::imag() const {
   return this->i_;
}

Complex Complex::operator+(Complex const& c2) const {
   return Complex(this->real() + c2.real(), this->imag() + c2.imag());
}

Complex Complex::operator-(Complex const& c2) const {
   return Complex(this->real() - c2.real(), this->imag() - c2.imag());
}

Complex Complex::operator*(Complex const& c2) const {
   return Complex(this->real() * c2.real() - this->imag() * c2.imag(),
                  this->real() * c2.imag() + this->imag() * c2.real());
}

Complex Complex::operator/(Complex const& c2) const {
   Complex conj = Complex(c2.real(), -c2.imag());
   Complex tmp = *this * conj;
   double norm = (c2 * conj).real();
   return Complex(tmp.real() / norm, tmp.imag() / norm);
}

bool Complex::operator==(Complex const& c2) const {
   Complex difference = *this - c2;
   return (std::abs(difference.real()) < Complex::EPSILON)
       && (std::abs(difference.imag()) < Complex::EPSILON);
}

bool Complex::operator!=(Complex const& c2) const {
   return !(*this == c2);
}

std::ostream& operator<<(std::ostream& out, const Complex& c) {
  if (c.imag() == 0) {
    out << c.real();
  }
  else if (c.real() == 0) {
    if (c.imag() == 1) {
      out << 'i';
    }
    else if (c.imag() == -1) {
      out << "-i";
    }
    else {
      out << c.imag() << 'i';
    }
  }
  else {
    if (c.imag() == 1) {
      out << c.real() << " + i";
    }
    else if (c.imag() == -1) {
      out << c.real() << " - i";
    }
    else if (c.imag() > 0) {
      out << c.real() << " + " << c.imag() << 'i';
    }
    else {
      out << c.real() << " - " << std::abs(c.imag()) << 'i';
    }
  }
  return out;
}
