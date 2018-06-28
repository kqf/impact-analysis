from cmath import pi, sqrt, exp
k_norm = 0.389379338


def amplitude(t):
    a1, a2, b1, b2, b3, b4, a_s, rho = [
        0.6335428499086242,
        1.6020821768950455,
        0.5841414884915991,
        1.0324336167108956,
        5.067007084850212,
        1.0030468301666557,
        1.10600e+02,
        0.1
    ]
    a_s = a_s / (sqrt(pi * k_norm) * 4)
    alpha = (1 - 1j * rho) * (a_s + a2)

    ampl = (
        1j * alpha * (
            exp(-0.5 * alpha * b1 * t) * a1 +
            exp(-0.5 * alpha * b2 * t) * (1 - a1)
        ) -
        1j * exp(-0.5 * b3 * t) * a2 - a2 * rho / ((1 + t / b4)**4))
    return ampl


def diff_cs(t):
    return abs(amplitude(t)) ** 2


def main():
    print diff_cs(0.010303)


if __name__ == '__main__':
    main()
