#include <cstdlib>
#include <cassert>
#include <iterator>
#include <vector>
#include <exception>
#include <memory>

using namespace std;

class OverflowExc : public ::std::exception {};

template <typename T>
class VecOb {
public:
	typedef vector<T> vec_t;
	typedef typename vec_t::const_reference const_reference;
	typedef typename vec_t::value_type value_type;
	typedef typename vec_t::difference_type difference_type;
	typedef typename vec_t::pointer pointer;
	typedef typename vec_t::reference reference;
	typedef typename vec_t::iterator::iterator_category iterator_category;
	vec_t data;

	VecOb() :
		data(1, T()) {}

	VecOb(typename vec_t::size_type maxWanted) :
		data(maxWanted + 1, T()) {}

	VecOb(typename vec_t::size_type maxWanted, const T &a) :
		data(maxWanted + 1, a) {}

	typename vec_t::size_type OnePast() const {
		return data.size();
	}

	const T & operator[](size_t a) const {
		assert(a >= 1);
		assert(a < OnePast());
		return data.at(a);
	}

	T & operator[](size_t a) {
		assert(a >= 1);
		assert(a < OnePast());
		return data.at(a);
	}

	typename vec_t::size_type      size() const  { return OnePast() - 1; } /* Size of [x, y) is 'y-x'. Here: [1, OnePast) */
	typename vec_t::const_iterator begin() const { return ++data.begin(); }
	typename vec_t::const_iterator end() const { return data.end(); }
	void                           push_back(const T& val) { data.push_back(val); }
};

namespace OneD {
	const int dimRowImage = 5;

	VecOb<int> MergePrefix(int i1, const VecOb<int> &j) {
		VecOb<int> w;

		w.push_back(i1);
		for (auto &i : j)
			w.push_back(i);

		return w;
	}

	VecOb<int> PasteInfix(const VecOb<int> &pre, int at, const VecOb<int> &post) {
		VecOb<int> w;

		for (auto &i : pre)
			w.push_back(i);
		w.push_back(at);
		for (auto &i : post)
			w.push_back(i);

		return w;
	}

	VecOb<int> CutSuffix(int d, const VecOb<int> &j) {
		VecOb<int> w;

		for (size_t i = d; i < j.OnePast(); i++)
			w.push_back(j[i]);

		return w;
	}

	/* [1, n_d] indexing */
	struct VoxelRef {
		bool undef;
		VecOb<int> i;

	private:
		VoxelRef(int d) : undef(true), i(d, 0) {}

	public:

		void Check() {
			if (undef) {
				for (auto &w : i) assert(w == 0);
			} else {
				for (auto &w : i) assert(w >= 1);
				/* FIXME: Unknown upper bound: for (auto &i : index) assert(i <= ???); */
			}
		}

		void SetFrom(const VoxelRef &v) {
			assert(i.size() == v.i.size());

			undef = v.undef;
			for (size_t i = 1; i < v.i.OnePast(); i++)
				Set(i, v.i[i]);

			Check();
		}

		void Set(int d, int newI) {
			i[d] = newI;
			Check();
		}

		static VoxelRef DOf(const VoxelRef &vOther, int dToSet, int newI) {
			VoxelRef w(vOther.i.size());
			w.SetFrom(vOther);
			w.Set(dToSet, newI);
			w.Check();
			return w;
		}

		static VoxelRef MakeUndef() {
			return VoxelRef(0);
		}

		static VoxelRef MakeMergingPrefix(int i1, const VecOb<int> &j) {
			VoxelRef ret(1 + j.size());
			VecOb<int> corVec = MergePrefix(i1, j);

			assert(ret.i.size() == corVec.size());

			ret.undef = false;

			for (size_t w = 1; w < corVec.OnePast(); w++)
				ret.Set(w, corVec[w]);

			return ret;
		}

		static VoxelRef MakeFromVec(const VecOb<int> &v) {
			VoxelRef ret(v.size());

			ret.undef = false;

			for (size_t w = 1; w < v.OnePast(); w++)
				ret.Set(w, v[w]);

			return ret;
		}
	};

	struct N {
		VecOb<int> dims;
		VecOb<int> accDims;

		N(const VecOb<int> &dms) :
			dims(dms)
		{
			assert(dms.size() >= 1);

			accDims = VecOb<int>(dms.size());
			accDims[1] = 1;
			for (size_t i = 2; i < accDims.OnePast(); i++)
				accDims[i] = accDims[i - 1] * dims[i - 1];
		}

		int & operator[](const size_t i) { return dims[i]; }
	};

	/* dims: [1, n_d] */
	struct I {
		N dims;
		vector<int> data;

		I(const VecOb<int> &dms) :
			dims(dms)
		{
			/* Do the [0, OnePast) overallocation. '(OnePast aka ''dims[i] + 1'') - 0' */
			int total = 1;
			for (size_t i = 1; i < dims.dims.OnePast(); i++)
				total *= (dims[i] + 1) - 0;
			data = vector<int>(total, 0);
		}

		void Check() {
			assert(dims.dims.size() == dims.accDims.size());
		}

		int GetAtRaw(const VecOb<int> &corVec) {
			assert(dims.dims.size() == corVec.size());

			int pos = 0;
			for (size_t i = 1; i < corVec.OnePast(); i++)
				pos += dims.accDims[i] * corVec[i];

			return data.at(pos);
		}

		int At(const VoxelRef &v) {
			return GetAtRaw(v.i);
		}

		/* d: varying coordinate; r: rest */
		vector<VoxelRef> GetRow(int d, const VoxelRef &r) {
			vector<VoxelRef> vr;

			for (int i = 1; i <= dims[d]; i++)
				vr.push_back(VoxelRef::DOf(r, d, i));

			return vr;
		}
	};

	struct F {
		N dims;
		vector<VoxelRef> data;

		F(const VecOb<int> &dms) :
			dims(dms)
		{
			int total = 1;
			for (size_t i = 1; i < dims.dims.OnePast(); i++) total *= (dims[i] + 1) - 0;
			data = vector<VoxelRef>(total, VoxelRef::MakeUndef());
		}

		size_t GetIdxOf(const VoxelRef &at) const {
			assert(dims.dims.size() == at.i.size());

			int pos = 0;
			for (size_t i = 1; i < at.i.OnePast(); i++)
				pos += dims.accDims[i] * at.i[i];

			return pos;
		}

		void Set(const VoxelRef &at, const VoxelRef &to) {
			data.at(GetIdxOf(at)) = to;
		}

		VoxelRef & operator[](const VoxelRef &a) {
			return data.at(GetIdxOf(a));
		}
	};

	namespace VecCount {

		void ReinitLower(int d, VecOb<int> *vInOut) {
			for (size_t i = d - 1; i >= 1; i--)
				(*vInOut)[i] = 1;
		}

		void Inc_(int d, const VecOb<int> &limits, VecOb<int> *vInOut) {
			if (d == vInOut->OnePast())
				throw OverflowExc();

			(*vInOut)[d] += 1;
			if ((*vInOut)[d] > limits[d])
				Inc_(d + 1, limits, vInOut);
			else
				ReinitLower(d, vInOut);
		}

		bool Inc(const VecOb<int> &limits, VecOb<int> *vInOut) {
			assert(limits.size() == vInOut->size());
			bool isOverflow = false;
			try { Inc_(1, limits, vInOut); } catch (const OverflowExc &) { isOverflow = true; }
			return isOverflow;
		}

		VecOb<int> MakeInitial(const VecOb<int> &limits) {
			return VecOb<int>(limits.size(), 1);
		}

	};
};

using namespace OneD;

struct Maurer {
	N varN;
	I varI;
	F varF;

private:
	Maurer(const VecOb<int> &dms) :
		varN(dms),
		varI(dms),
		varF(dms) {}

public:
	static Maurer * MakeUniform(int d, int n) {
		VecOb<int> dms;
		for (int i = 1; i <= d; i++) dms.push_back(n);
		return new Maurer(dms);
	}

	int EucDist(const VoxelRef &a, const VoxelRef &b) {
		assert(0);
		return 0;
	}

	bool RemoveFT(const VoxelRef &glm1, const VoxelRef &gl, const VoxelRef &fi, void *RowD) {
		assert(0);
		return false;
	}

	void VoronoiFT(int d, const VecOb<int> &is, const VecOb<int> &js) {
		int l = 0;
		VecOb<VoxelRef> gs(varN.dims[d], VoxelRef::MakeUndef());

		for (int i = 1; i <= varN.dims[d]; i++) {
			VoxelRef xi = VoxelRef::MakeFromVec(PasteInfix(is, i, js));
			VoxelRef fi = varF[xi];

			if (! fi.undef) {
				if (l < 2) {
					gs[++l] = fi;
				} else {
					while (l >= 2 && RemoveFT(gs[l - 1], gs[l], fi, nullptr))
						--l;
					gs[++l] = fi;
				}
			}
		}

		int ns = l;
		l = 1;

		if (ns == 0)
			return;

		for (int i = 1; i <= varN.dims[d]; i++) {
			VoxelRef xi = VoxelRef::MakeFromVec(PasteInfix(is, i, js));
			while (l < ns && EucDist(xi, gs[l]) > EucDist(xi, gs[l + 1]))
				++l;
			varF[xi] = gs[l];
		}
	}

	void ComputeFT(int d, const VecOb<int> &j) {
		if (d == 1) {
			assert(d + j.size() == varN.dims.size());

			for (int i1 = 1; i1 <= varN.dims[1]; i1++) {
				VoxelRef w = VoxelRef::MakeMergingPrefix(i1, j);
				VoxelRef u = VoxelRef::MakeUndef();
				if (varI.At(w) == 1)
					varF[w] = w;
				else
					varF[w] = u;
			}
		} else {
			for (int id = 1; id <= varN.dims[d]; id++)
				ComputeFT(d - 1, MergePrefix(id, j));
		}

		/* FIXME: Boundary case. i1...id-1 where d == 1
		Not sure what the intention is for d == 1.
		After the initial part computes F_(d-1), the routine should compute F_d.
		Maybe VoronoiFT can be dummy-called (Empty 'i' part). */
		const int lastNested = d - 1;
		VecOb<int> limits;
		copy_n(varN.dims.begin(), lastNested, back_inserter(limits));

		VecOb<int> count = VecCount::MakeInitial(limits);

		do {
			VoronoiFT(d, count, CutSuffix(d + 1, j)); /* FIXME: Pass the correct 'j' */
		} while (!VecCount::Inc(limits, &count));
	}

	void Start() {
		ComputeFT(varN.dims.size(), VecOb<int>(0));
	}
};

int main(int argc, char **argv) {
	shared_ptr<Maurer> m(Maurer::MakeUniform(1, 5));
	m->Start();
	return EXIT_SUCCESS;
}
