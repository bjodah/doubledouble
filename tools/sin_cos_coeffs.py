import mpmath


mpmath.mp.dps = 100


def split_to_doubledouble(x):
    upper = float(x)
    lower = float(x - mpmath.mpf(upper))
    return upper, lower


def print_doubledouble_array(name, values):
    print(f"static const std::array<DoubleDouble, {len(values)}> {name}{{{{")
    for k, value in enumerate(values, start=1):
        upper, lower = split_to_doubledouble(value)
        comma = "," if k < len(values) else ""
        print(f"    {{{upper: .17e}, {lower: .17e}}}{comma}  // k={k}")
    print("}};")


def sin_values(denom, nterms):
    return [mpmath.sin(k * mpmath.pi / denom) for k in range(1, nterms + 1)]


def cos_values(denom, nterms):
    return [mpmath.cos(k * mpmath.pi / denom) for k in range(1, nterms + 1)]


def print_table_group(denom, nterms, suffix):
    print(f"// Table of sin(k*pi/{denom}) and cos(k*pi/{denom}), k=1..{nterms}.")
    print_doubledouble_array(f"sin_table_{suffix}", sin_values(denom, nterms))
    print()
    print_doubledouble_array(f"cos_table_{suffix}", cos_values(denom, nterms))


print("// qd-compatible check tables. These should reproduce dd_real.cpp's pi/16 table.")
print_table_group(16, 4, "pi_16")
print()
print("// Proposed implementation tables for doubledouble trig range reduction.")
print_table_group(24, 6, "pi_24")
print()
upper, lower = split_to_doubledouble(mpmath.pi / 24)
print(f"inline const DoubleDouble dd_pi_24{{{upper: .17e}, {lower: .17e}}};")
