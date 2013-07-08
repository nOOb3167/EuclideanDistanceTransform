#include <cstdlib>
#include <cassert>
#include <iterator>
#include <vector>

using namespace std;

template <typename T>
class VecOb {
public:
	typedef vector<T> vec_t;
	vec_t data;

	VecOb() :
		data(1, T()) {}

	VecOb(vec_t::size_type maxWanted) :
		data(maxWanted + 1, T()) {}

	VecOb(vec_t::size_type maxWanted, const T &a) :
		data(maxWanted + 1, a) {}

	vec_t::size_type OnePast() const {
		return data.size() + 1;
	}

	T & operator[](size_t a) const {
		assert(a >= 1);
		assert(a < OnePast());
		return data.at(a);
	}

	vec_t::size_type      size() const  { return OnePast() - 1; } /* Size of [x, y) is 'y-x'. Here: [1, OnePast) */
	vec_t::const_iterator begin() const { return data.begin() ++; }
	vec_t::const_iterator end() const   { return data.end(); }
	void                  push_back(const T& val) { data.push_back(val); }
};

namespace OneD {
	const int dimRowImage = 5;

	VecOb<int> MergePrefix(int i1, const VecOb<int> &j) {
		VecOb<int> w;

		w.push_back(i1);
		for (auto &i : j) w.push_back(i);

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
				for (auto &i : index) assert(i == 0);
			} else {
				for (auto &i : index) assert(i >= 1);
				/* FIXME: Unknown upper bound: for (auto &i : index) assert(i <= ???); */
			}
		}

		void SetFrom(const VoxelRef &v) {
			assert(i.size() == v.i.size());

			undef = v.undef;
			for (size_t i = 1; i < v.i.OnePast(); i++)
				Set(i, v.i[i]);
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
			VoxelRef w(1 + j.size());
			VecOb<int> corVec = MergePrefix(i1, j);

			assert(w.i.size() == corVec.size());

			for (size_t i = 1; i < w.i.OnePast(); i++)
				w.Set(i, corVec[i]);

			return w;
		}
	};

	struct F {
		/* FIXME: Empty */
	};

	struct N {
		VecOb<int> dims;
		VecOb<int> accDims;

		N(const VecOb<int> &dms) :
			dims(dms)
		{
			assert(dms.size() >= 1);

			VecOb<int> wAccDims(dms.size());
			wAccDims[1] = 1;
			for (int i = 2; i < wAccDims.OnePast(); i++) wAccDims[i] = wAccDims[i - 1] * wDims[i - i];

			accDims = wAccDims;
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
			vector<int> wData;
			{
				int total = 1;
				for (size_t i = 1; i < dims.dims.OnePast(); i++)
					total *= (dims[i] + 1) - 0;
				wData = vector<int>(total, 0);
			}

			data = wData;
		}

		void Check() {
			assert(dims.dims.size() == dims.accDims.size());
		}

		int GetAtRaw(const VecOb<int> &corVec) {
			assert(dims.dims.size() == corVec.size());

			int pos = 0;
			for (int i = 1; i < corVec.OnePast(); i++)
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
};

using namespace OneD;

struct Maurer {
	N varN;
	I varI;
	F varF;

	void ComputeFT(int d, const VecOb<int> &j) {
		if (d == 1) {
			assert(d + j.size() == varN.dims.size());

			for (int i1 = 1; i1 < varN.dims.OnePast(); i1++) {
				VoxelRef w = VoxelRef::MakeMergingPrefix(i1, j);
				VoxelRef u = VoxelRef::MakeUndef();
				if (varI.At(w) == 1)
					varF.Set(w, w);
				else
					varF.Set(w, u);
			}
		} else {
			for (int id = 1; id <= varN.At(d); id++)
				ComputeFT(d - 1, MergePrefix(id, j));
		}

		/* FIXME: Boundary case. i1...id-1 where d == 1 */
		int numNested;
		if (d == 1) numNested = 1;
		else        numNested = d - 1;

		vector<int> tmpId(numNested, 1);
	}
};

int main(int argc, char **argv) {
	return EXIT_SUCCESS;
}
