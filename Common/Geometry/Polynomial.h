#pragma once
#include <vector>
#include "FEMlibrary.h"
using namespace std;

namespace FemLibrary {

	class FEMLIBRARY_API Polynomial
	{
	public:


		template<typename T>
		static vector<T> FindRootsQuadratic(const vector<T>& poly)
		{
			vector<T> res;
			assert(poly.size() == 2); // x^2+bx+c=0
			auto c = poly[0];
			auto b = poly[1];
			auto d = b*b - 4 * c;
			if (d < 0) return res;

			auto sqrtD = sqrt(d);
			res.resize(2);
			res[0] = (-b - sqrtD) / 2;
			res[1] = (-b + sqrtD) / 2;
			if (res[0] == res[1]) res.resize(1); // явное сравнение
			return res;
		}

#define M_PI  (3.14159265358979323846)
#define M_2PI (2.*M_PI)

		// https://en.wikipedia.org/wiki/Cubic_function
		// http://algolist.manual.ru/maths/findroot/cubic.php
		template<typename T>
		static vector<T> FindRootsCubic(const vector<T>& poly)
		{
			vector<T> res;
			assert(poly.size() == 3); // x^3 + a x^2 + b x + c = 0
			auto c = poly[0];
			auto b = poly[1];
			auto a = poly[2];

			auto a2 = a*a;
			auto q = (a2 - 3.*b) / 9.;
			auto r = (a*(2.*a2 - 9.*b) + 27.*c) / 54.;
			auto r2 = r*r;
			auto q3 = q*q*q;
			auto s = q3 - r2;
			a /= 3.;
			if (s >= 0)
			{
				auto sqrtq3 = sqrt(q3);
				if (sqrtq3 == 0)
					res.emplace_back(-a);
				else
				{
					res.resize(3);
					auto t = acos(r / sqrtq3); // todo: сделать сравнение по r / sqrt(q3)
					q = -2.*sqrt(q);
					res[0] = q*cos(t / 3.) - a;
					res[1] = q*cos((t + M_2PI) / 3.) - a;
					res[2] = q*cos((t - M_2PI) / 3.) - a;
					assert(res[0] != res[1]);
					if (res[0] == res[2] || res[1] == res[2]) res.resize(2); // явное сравнение
				}
			}
			else
			{
				auto signR = r >= 0. ? 1.0 : -1.0;
				auto aa = -signR * pow(signR * r + sqrt(-s), 1. / 3.); // A = -sign(R) [|R|+sqrt(R2-Q3)] ^ 1/3
				auto bb = aa != 0. ? q / aa : 0.;
				res.emplace_back(aa + bb - a);
			}
			return res;
		}

		// https://en.wikipedia.org/wiki/Cubic_function
		// http://algolist.manual.ru/maths/findroot/cubic.php
		template<typename T>
		static vector<T> FindRootsQuartic(const vector<T>& poly)
		{
			vector<T> res;
			assert(poly.size() == 4); //  x^4 + b x^3 + c x^2 + d x + e = 0
			auto e = poly[0];
			auto d = poly[1];
			auto c = poly[2];
			auto b = poly[3];

			auto b2 = b * b;
			auto c2 = c * c;
			auto p = c - (3. / 8.) * b2;
			auto q = b * (b2 - 4 * c) / 8 + d;
			//auto r = (b2 * (16 * c - 3 * b2) - 64 * b * d) / 256 + e;
			auto delta0 = c2 - 3 * b*d + 12 * e;
			auto delta1 = c(2 * c2 - 9 * b*d - 72 * c*e) + 27 * (b2*e + d*d);
			auto delta12 = delta1*delta1;
			auto delta03 = delta0*delta0*delta0;
			auto delta = delta12 - 4 * delta03;
			auto discriminant = -delta / 27;

			if (discriminant > 0)
			{
				auto bigD = b2*(16 * c - 3 * b2) - 16 * (b*d + c2) + 64 * e;
				if (p > 0 || bigD > 0) return res; // no real roots

				auto t = acos(delta1 / 2 / sqrt(delta03)) / 3;
				auto s = sqrt(2. / 3.*(sqrt(delta0)*cos(t) - p)) / 2;
				res.resize(4);

				auto root = -4 * s * s - 2 * p + q / s;
				assert(root > 0);
				root = sqrt(root) / 2;
				res[0] = -s - root;
				res[1] = -s + root; // todo: sort 1 and 2
				res[2] = s - root;
				res[3] = s + root;
			}
			else
			{
				auto Q3 = (delta1 + sqrt(delta)) / 2;
				auto signQ = Q3 < 0 ? -1 : 1;
				auto Q = signQ * pow(signQ * Q3, 1. / 3.);
			}

			for (auto& root : res)
				root -= b / 4;
			return res;
		}

		static constexpr double xeps = 1e-6;
		static constexpr double yeps = 1e-10;
		static constexpr int newtons_method_max_iterations = 100;

		static bool NewtonsMethod(const vector<double>& poly, const vector<double>& derivative, double a, double b, double& root);
		static double EvaluatePolynomial(const vector<double>& poly, double x);
		static vector<double> SolveQuinticInRange(vector<double>& poly, double a, double b, bool max_root_only = false);
		static vector<double> FindRoots(const vector<double>& poly);
		static double CheckPolynomialSolution(int n);
		static double CheckQuinticSolution();
	};

}