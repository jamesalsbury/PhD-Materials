#include <Module.h> // include JAGS module base class
2 #include <distributions/DBern.h> // include Bernoulli distribution class
3
4 namespace Bernoulli { // start defining the module namespace
5
6 // Module class
7 class BERNModule : public Module {
8 public:
9 BERNModule(); // constructor
10 ~BERNModule(); // destructor
11 };
12
13 // Constructor function
14 BERNModule::BERNModule() : Module("Bernoulli") {
15 insert(new DBern); // inherited function to load objects into JAGS
16 }
17
18 // Destructor function
19 BERNModule::~BERNModule() {
20 std::vector<Distribution*> const &dvec = distributions();
21 for (unsigned int i = 0; i < dvec.size(); ++i) {
22 delete dvec[i]; // delete all instantiated distribution objects
23 }
24 }
25
26 } // end namespace definition
27
28 Bernoulli::BERNModule _Bernoulli_module;