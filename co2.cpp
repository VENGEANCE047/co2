#include <bits/stdc++.h>
using namespace std;

double mi_(double wi)
{
    double mi = 0.37464 + 1.54226 * wi - 0.26992 * wi * wi;
    return mi;
}
double a_(double R, double T, double Tc_co2, double Pc_co2, double mi)
{
    double a = (0.457236 * R * R * Tc_co2 * Tc_co2 / Pc_co2) * pow((1 + mi * (1 - pow(T / Tc_co2, 0.5))), 2);
    return a;
}
double b_(double R, double Tc_co2, double Pc_co2)
{
    double b = 0.077796 * R * Tc_co2 / Pc_co2;
    return b;
}
double A_(double a, double P, double R, double T)
{
    double A = a * P / (R * R * T * T);
    return A;
}
double B_(double b, double P, double R, double T)
{
    double B = b * P / (R * T);
    return B;
}
double f_(double Z_old, double A, double B)
{
    double f = Z_old * Z_old * Z_old - (1 - B) * Z_old * Z_old + (A - 2 * B - 3 * B * B) * Z_old - (A * B - B * B - B * B * B);
    return f;
}
double df_(double Z_old, double A, double B)
{
    double df = 3 * Z_old * Z_old - 2 * Z_old * (1 - B) + (A - 2 * B - 3 * B * B);
    return df;
}
double Z_(double A, double B)
{
    double Z_old;
    int c1 = 0;
    double Z = Z_old;
    do
    {
        c1++;
        Z_old = Z;
        Z = Z_old - f_(Z_old, A, B) / df_(Z_old, A, B);
    } while (abs(Z - Z_old) > 0.000001);
    return Z;
}

double phii_(double Z, double B, double A, double del2, double del1, double a)
{
    double phii = pow(M_E, ((Z - 1) - log(Z - B) - (A / (B * (del2 - del1)) * (log((Z + del2 * B) / (Z + del1 * B))))));
    return phii;
}
double lamna_(double T, double P)
{
    double lamna = (1.6790636 * pow(10, -4) * T) + (40.838951 / T) - 0.0652869 - (3.9266518 * pow(10, -2) * P / T) + (2.1157167 * pow(10, -2) * P / (630 - T)) + (6.5486487 * pow(10, -6) * T * log(P));
    return lamna;
}
double epnacl_(double T, double P)
{
    double epnacl = (2.8274958 * pow(10, -5) * T) - (1.144624 * pow(10, -2)) + ((1.3980876 * pow(10, -2) * P) / T) - (1.4349005 * pow(10, -2) * P) / (630 - T);
    return epnacl;
}
double gammai_(double mc, double lamna, double epnacl)
{
    double gammai = exp((2 * mc * lamna) + (2 * mc * mc * epnacl));
    return gammai;
}
double delB_(double T)
{
    double delB = -5.279063 + 6.187967 * pow((1000.00 / T), 0.5);
    return delB;
}
double Ps_(double Pc_h2o, double Tc_h2o, double T)
{
    double Ps = Pc_h2o * pow(M_E, (Tc_h2o / T) * (1.8440825 * pow((1 - T / Tc_h2o), 1.5) - 7.8595178 * (1 - T / Tc_h2o) - 11.786649 * pow((1 - T / Tc_h2o), 3) + 22.680741 * pow((1 - T / Tc_h2o), 3.5) - 15.9618719 * pow((1 - T / Tc_h2o), 4) + 1.8012250 * pow((1 - T / Tc_h2o), 7.5)));
    return Ps;
}
double V0_(double theta)
{
    double V0 = (1 + (18.1597 * 1e-3 * theta)) / (0.9998 + (18.2249 * 1e-3 * theta) - (7.9222 * 1e-6 * theta * theta) - (55.4485 * 1e-9 * pow(theta, 3)) + (149.7562 * 1e-12 * pow(theta, 4)) - (393.2952 * 1e-15 * pow(theta, 5)));
    return V0;
}
double A1_(double theta)
{
    double A1 = 3.2891 - 2.391 * 0.001 * theta + 2.8446 * 0.0001 * theta * theta - 2.82 * pow(10, -6) * theta * theta * theta + 8.447 * pow(10, -9) * pow(theta, 4);
    return A1;
}
double A2_(double theta)
{
    double A2 = 6.245 * pow(10, -5) - 3.913 * pow(10, -6) * theta - 3.499 * pow(10, -8) * pow(theta, 2) + 7.942 * pow(10, -10) * pow(theta, 3) - 3.299 * pow(10, -12) * pow(theta, 4);
    return A2;
}
double B1_(double theta)
{
    double B1 = 19654.32 + 147.037 * theta - 2.2155 * pow(theta, 2) + 1.0478 * 0.01 * pow(theta, 3) - 2.2789 * pow(10, -5) * pow(theta, 4);
    return B1;
}
double rho_(double V0, double A1, double A2, double B1, double P)
{
    double rho = (B1 + (A1 * P) + (A2 * P * P)) / (B1 * V0 + (A1 * P * V0) + (A2 * P * P * V0) - (V0 * P));
    return rho;
}
double f0h2o_(double Ps, double P, double rho, double T, double R)
{
    double f0h2o = Ps * pow(M_E, 18.0152 * (P - Ps) / (rho * R * T));
    return f0h2o;
}
double hi_(double T, double eta, double f0h2o, double rho0h2o, double mwh2o, double R, double delB)
{
    double hi = ((1 - eta) * log(f0h2o) + (eta * log(R * T * rho0h2o / mwh2o)) + (2 * rho0h2o * delB));
    return exp(hi);
}
double Ki_(double hi, double gammai, double P, double phii)
{
    double Ki = (hi * gammai) / (P * phii);
    return Ki;
}
double K0h2o_(double theta)
{
    double term = (-2.209) + (3.097 * (1e-2) * theta) - (1.098 * (1e-4) * theta * theta) + (2.048 * (1e-7) * theta * theta * theta);
    double K0h2o = pow(10, term);
    return K0h2o;
}

double Kh2o_(double K0h2o, double f0h2o, double P, double R, double T)
{
    double Kh2o = (K0h2o / (f0h2o * P)) * pow(M_E, (((P - 1) * 18.18) / (R * T)));
    return Kh2o;
}
double yh2o_(double Ki, double Kh2o)
{
    double yh2o = (1.0 - 1.0 / Ki) / (1.0 / Kh2o - 1.0 / Ki);
    return yh2o;
}

double yin_(double yi, double yh2o)
{
    double yin = yi / (1 + yh2o);
    return yin;
}
double xi_(double yin, double Ki)
{
    double xi = yin / Ki;
    return xi;
}

int main()
{
    double yi = 1;
    double P, phii, hi, Ki, yh2o, yin, xi, b, a, A, B, mc, gammai, lamna, epnacl, delB;
    double del2 = 1 - pow(2, 0.5);
    double del1 = 1 + pow(2, 0.5);
    double eta = -0.114535;
    double f0h2o;
    double rho0h2o;
    double mwh2o = 18.015;
    double R = 83.14;
    double Kh2o;
    double K0h2o;
    double T;
    double theta;
    double Tc_h2o = 374 + 273.15;
    double Pc_h2o = 220.64;
    double Tc_co2 = 30.9780 + 273.15;
    double Pc_co2 = 73.8;
    double mi;
    double wi = 0.224;
    double Ps;
    double rho;
    double A1, A2, B1, V0;
    double Z_old = 1;
    double Z = 2;
    double f, df;

    ifstream inputFile("input.txt");
    if (!inputFile)
    {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    ofstream outputFile("Output.txt");
    while (inputFile >> T >> P >> mc)
    {
        theta = T - 273.15;
        mi = mi_(0.224);
        a = a_(R, T, Tc_co2, Pc_co2, mi);
        b = b_(R, Tc_co2, Pc_co2);
        B = B_(b, P, R, T);
        A = A_(a, P, R, T);
        f = f_(Z_old, A, B);
        df = df_(Z_old, A, B);
        Z = Z_(A, B);
        phii = phii_(Z, B, A, del2, del1, a);
        lamna = lamna_(T, P);
        epnacl = epnacl_(T, P);
        gammai = gammai_(mc, lamna, epnacl);
        delB = delB_(T);
        Ps = Ps_(Pc_h2o, Tc_h2o, T);
        V0 = V0_(theta);
        A1 = A1_(theta);
        A2 = A2_(theta);
        B1 = B1_(theta);
        rho = rho_(V0, A1, A2, B1, P);
        f0h2o = f0h2o_(Ps, P, rho, T, R);
        hi = hi_(T, eta, f0h2o, rho, mwh2o, R, delB);
        Ki = Ki_(3515, gammai, P, phii);
        K0h2o = K0h2o_(theta);
        Kh2o = Kh2o_(K0h2o, f0h2o, P, R, T);
        yh2o = yh2o_(Ki, Kh2o);
        yin = yin_(yi, yh2o);
        xi = xi_(yin, Ki);
        outputFile << "Fugacity Coefficient for CO2 : " << phii << endl;
        outputFile << "Henry's coefficient for CO2 : " << hi << endl;
        outputFile << "Equilbrium constant for H2O : " << Ki << endl;
        outputFile << "Vapour Phase Mole Fraction of H2O : " << yin << endl;
        outputFile << "Liquid Phase Mole Fraction of CO2 : " << xi << endl;
        outputFile << endl;
    }
}