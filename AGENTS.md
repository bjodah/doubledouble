This is a stand-alone C++ library for double-double arithmetics.

## Testing
To compile and run tests:
```console
$ cd tests
$ make
$ ./test_doubledouble
```

When designing tests, it is important to derive high quality reference data for numerical values.
One simple way is to use e.g. `mpmath` in Python:
```console
$ python - <<'PYEOF'
> import mpmath
with mpmath.workdps(70):
    atan1 = mpmath.atan(1)
    four_atan1 = 4*atan1
    four_atan1_f64_high = float(four_atan1)
    four_atan1_f64_low  = float(four_atan1 - four_atan1_f64_high)
    print(f"{four_atan1_f64_high:23.17e}", end=", ")
    print(f"{four_atan1_f64_low:23.17e}")
PYEOF
3.14159265358979312e+00, 1.22464679914735321e-16
```
For benchmarking the implementation we can compare with boost...
```console
$ cat boost_dd.cpp
#include <cstdio>
#include <boost/version.hpp>
#if BOOST_VERSION < 108900
#error "Need boost 1.89 or newer for doubledouble"
#endif
#include <boost/multiprecision/cpp_double_fp.hpp>
using bmpdd = boost::multiprecision::cpp_double_double;

int main() {
    bmpdd one {1};
    bmpdd atan1 = atan(one);
    bmpdd four_atan1 = bmpdd{4} * atan1;
    double four_atan1_f64_high = double{four_atan1};
    double four_atan1_f64_low = double{four_atan1 - four_atan1_f64_high};
    std::printf("%23.17e, %23.17e\n", four_atan1_f64_high, four_atan1_f64_low);
    return 0;
}
$ c++ -isystem /opt-4/boost-1.91.0-release/include -O2 -Wall -Wextra -o boost_dd boost_dd.cpp && ./boost_dd
3.14159265358979312e+00, 1.22464679914735296e-16
```
...and libqd:
```console
$ grep -E '^struct .* dd_real' /src-2/qd-2.3.24/include/qd/dd_real.h
struct QD_API dd_real {
$ grep -E 'dd_real::dd_real' /src-2/qd-2.3.24/src/dd_real.cpp | head -1
dd_real::dd_real(const char *s) {
$ 
$ cat libqd_dd.cpp
#include <cstdio>
#include <qd/qd_real.h>

int main() {
    dd_real one {1};
    dd_real atan1 = atan(one);
    dd_real four_atan1 = dd_real{4} * atan1;
    double four_atan1_f64_high = to_double(four_atan1);
    double four_atan1_f64_low = to_double(four_atan1 - four_atan1_f64_high);
    std::printf("%23.17e, %23.17e\n", four_atan1_f64_high, four_atan1_f64_low);
    return 0;
}
$ c++ -isystem /opt-2/qd-2.3.24-release/include -O2 -Wall -Wextra -o libqd_dd libqd_dd.cpp -L/opt-2/qd-2.3.24-release/lib -Wl,-rpath,/opt-2/qd-2.3.24-release/lib -lqd && ./libqd_dd
3.14159265358979312e+00, 1.22464679914735321e-16
```
