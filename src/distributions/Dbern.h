#ifndef DBERN_H_
2 #define DBERN_H_
3 #include <distribution/ScalarDist.h> // JAGS scalar distribution base
class
4
5 namespace Bernoulli {
6
7 class DBern : public ScalarDist // scalar distribution class
8 {
9 public:
10 DBern(); // constructor
11 double logDensity(double x, PDFType type,
12 std::vector<double const *> const &parameters,
13 double const *lower, double const *upper) const;
14 double randomSample(std::vector<double const *> const &parameters,
15 double const *lower, double const *upper,
16 RNG *rng) const;
17 double typicalValue(std::vector<double const *> const &parameters,
18 double const *lower, double const *upper) const;
19 /** Checks that pi lies in the open interval (0,1) */
20 bool checkParameterValue(std::vector<double const *> const &
parameters) const;
21 /** Bernoulli distribution is discrete valued */
22 bool isDiscreteValued(std::vector<bool> const &mask) const;
23 };
24
25 }
26 #endif /* DBERN_H_ */
