#include "includes.h"
#include "Polynomial.h"
#include "P51 root finder.h"

namespace FemLibrary {

	bool Polynomial::NewtonsMethod(const vector<double>& poly, const vector<double>& derivative, double a, double b, double& root)
	{
		auto x0 = (a + b) / 2;
		for (auto iteration = 0; iteration < newtons_method_max_iterations; ++iteration)
		{
			auto fx0 = EvaluatePolynomial(poly, x0);
			auto x1 = x0 - fx0 / EvaluatePolynomial(derivative, x0);
			while (x1 < a || x1 > b) x1 = (x0 + x1) / 2;
			if (abs(x0 - x1) < xeps)
			{
				root = (x0 + x1) / 2;
				if (abs(EvaluatePolynomial(poly, root)) < yeps)
					return true;
			}
			x0 = x1;
		}
		return false;

	}

	double Polynomial::EvaluatePolynomial(const vector<double>& poly, double x)
	{
		auto value = poly[0];
		for (size_t i = 1; i < poly.size(); ++i)
			value = poly[i] + x*value;
		return value;
	}

	vector<double> Polynomial::SolveQuinticInRange(vector<double>& poly, double a, double b, bool max_root_only)
	{	// старшие степени вперед

		assert(poly.size() == 6);

		for (auto i = 5; i >= 0; --i) // нормировка
			poly[i] /= poly[0];

		vector<double> derivative(5);
		for (auto i = 0; i < 5; i++)
			derivative[i] = (5 - i)*poly[i];

		auto extremums = FindRoots(derivative);
		extremums.emplace_back(a);
		sort(extremums.begin(), extremums.end());

		vector<double> solutions;
		auto x = b, x_value = EvaluatePolynomial(poly, b);
		auto x_old_value = x_value;

		for (ptrdiff_t i = extremums.size() - 1; i >= 0; --i)
		{
			if (extremums[i] >= x) continue;
			if (extremums[i] < a) break;

			auto new_value = EvaluatePolynomial(poly, extremums[i]);
			
			if (abs(x_value) < yeps && x_value*x_old_value >= 0 && new_value*x_value >= 0) { // abs(new_value) < yeps  // !!! явное сравнение
				solutions.emplace_back(x);
				solutions.emplace_back(x);
				if (max_root_only) return solutions;
			}
			
			if (x_value*new_value < 0)
			{
				double root;
				if (NewtonsMethod(poly, derivative, extremums[i], x, root))
				{
					solutions.emplace_back(root);
					if (max_root_only) return solutions;
				}
				else throw exception("SolveQuinticInRange: newton fails");
			}

			x_old_value = x_value;
			x_value = new_value;
			x = extremums[i];
		}
		if (abs(x_value) < yeps && x_value*x_old_value > 0) {
			solutions.emplace_back(x);
		}

		assert(solutions.size() <= 5);
		return solutions;
	}

	vector<double> Polynomial::FindRoots(const vector<double>& poly)
	{
		auto n = poly.size() - 1;
		vector<long double> kbuffer(n + 1);
		vector<long double> real_roots(n);
		vector<long double> imag_roots(n);

		for (size_t i = 0; i <= n; ++i)
			kbuffer[i] = poly[i];

		// Пятая степень в P51 глючит, поэтому используем только QuadCubicRoots
		//auto nroots = n < 5 ? QuadCubicRoots(kbuffer.data(), n, real_roots.data(), imag_roots.data()) 
		//	                : PFiftyOne		(kbuffer.data(), n, real_roots.data(), imag_roots.data());
		assert(n < 5);
		auto nroots = QuadCubicRoots(kbuffer.data(), static_cast<int>(n), real_roots.data(), imag_roots.data());

		// Этот эксепшн постоянно вылетает ибо условие неверное
		// if (nroots != n)
		// 	throw exception("can't solve polynomial");

		vector<double> res;
		for (int i = 0; i < nroots; ++i)
			if (imag_roots[i] == 0.0)
				res.emplace_back(real_roots[i]);
		return res;
	}

	double Polynomial::CheckPolynomialSolution(int n)
	{
		vector<long double> kbuffer(n + 1);
		vector<long double> real_roots(n);
		vector<long double> imag_roots(n, 0.0);
		vector<long double> real_roots2(n);
		vector<long double> imag_roots2(n);
		for (auto i = 0; i < n; ++i)
			real_roots[i] = 4 * static_cast<long double>(rand()) / RAND_MAX - 2; // [-2; 2]
		BuildPolyFromRoots(n, real_roots.data(), imag_roots.data(), kbuffer.data());

		// Пятая степень в P51 глючит, поэтому используем только QuadCubicRoots
		//auto nroots = n < 5 ? QuadCubicRoots(kbuffer.data(), n, real_roots2.data(), imag_roots2.data())
		//	                : PFiftyOne		(kbuffer.data(), n, real_roots2.data(), imag_roots2.data());
		assert(n < 5);
		auto nroots = QuadCubicRoots(kbuffer.data(), n, real_roots.data(), imag_roots.data());

		if (nroots == 0)
			return 666;// MAGIC
			//throw exception("can't solve polynomial");

		SortRoots(n, real_roots.data(), imag_roots.data(), real_roots2.data(), imag_roots2.data());

		auto mean_square_error = 0.0;
		for (auto i = 0; i < n; ++i)
		{
			auto delta = real_roots[i] - real_roots2[i];
			mean_square_error += delta*delta;
			delta = imag_roots[i] - imag_roots2[i];
			mean_square_error += delta*delta;
		}
		return mean_square_error / n / 2;
	}

	double Polynomial::CheckQuinticSolution()
	{
		auto n = 5;
		auto a = 0., b = 1.;
		vector<long double> kbuffer(n + 1);
		vector<long double> real_roots(n);
		vector<long double> imag_roots(n, 0.0);
		vector<double> roots;
		for (auto i = 0; i < n; ++i) {
			real_roots[i] = 4 * static_cast<long double>(rand()) / RAND_MAX - 2; // [-2; 2]
			if (a <= real_roots[i] && real_roots[i] <= b)
				roots.emplace_back(real_roots[i]);
		}
		BuildPolyFromRoots(n, real_roots.data(), imag_roots.data(), kbuffer.data());

		vector<double> poly(n + 1);
		transform(kbuffer.begin(), kbuffer.end(), poly.begin(), [](long double d) { return d; });
		auto roots2 = SolveQuinticInRange(poly, a, b);

		if (roots.size() != roots2.size())
			return 666; // MAGIC

		sort(roots.begin(), roots.end());
		sort(roots2.begin(), roots2.end());

		auto mean_square_error = 0.0;
		for (size_t i = 0; i < roots.size(); ++i)
		{
			auto delta = roots[i] - roots2[i];
			mean_square_error += delta*delta;
		}
		return mean_square_error / roots.size();
	}

}