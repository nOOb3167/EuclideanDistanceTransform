#include <cstdlib>
#include <cassert>
#include <iterator>
#include <vector>

using namespace std;

namespace OneD {
	const int dimRowImage = 5;

	/* [1, n_d] indexing */
	struct VoxelRef {
		bool undef;
		int i;
		VoxelRef() : undef(true), i(0) {}

		void Check() {
			if (undef) {
				assert(i == 0);
			} else {
				assert(i >= 1);
				assert(i <= 1);
			}
		}

		void SetFrom(const VoxelRef &v) {
			undef = v.undef;
			Set(1, v.i);
		}

		void Set(int d, int newI) {
			i = newI;
			Check();
		}

		int GetNumDims() const { return 1; }

		vector<int> GetCorVec() const {
			return vector<int>(1, i);
		}

		static VoxelRef DOf(const VoxelRef &vOther, int dToSet, int newI) {
			VoxelRef w;
			w.SetFrom(vOther);
			w.Set(dToSet, newI);
			w.Check();
			return w;
		}
	};

	struct F {
		VoxelRef data[dimRowImage];
	};

	/* dims: [1, n_d] */
	struct I {
		N dims;

		vector<int> data;

		I(const vector<int> dms) :
			dims(dms)
		{
			vector<int> wData;
			{
				int total = 1;
				for (int i = 1; i <= dims.GetNumDims(); i++) total *= dims.GetDim(i);
				wData = vector<int>(total, 0);
			}

			data    = wData;
		}

		void Check() {
			assert(dims.size() >= 1);
			assert(dims.size() == accDims.size());
		}

		int GetAtRaw(vector<int> corVec) {
			assert(GetNumDims() == corVec.size());

			int pos = 0;
			for (int i = 0; i < corVec.size(); i++) {
				pos += accDims[i] * (corVec[i] - 1);
			}

			return data.at(pos);
		}

		int At(const VoxelRef &v) {
			return GetAtRaw(v.GetCorVec());
		}

		/* d: varying coordinate; r: rest */
		vector<VoxelRef> GetRow(int d, const VoxelRef &r) {
			vector<VoxelRef> vr;

			for (int i = 1; i <= dims.GetDim(d); i++)
				vr.push_back(VoxelRef::DOf(r, d, i));

			return vr;
		}
	};

	struct N {
		vector<int> dims;
		vector<int> accDims;

		N(const vector<int> dms)
		{
			assert(dms.size() >= 1);

			vector<int> wDims(dms.size() + 1);
			for (int i = 1; i < wDims.size(); i++) wDims[i] = dms[i - 1];

			vector<int> wAccDims(dims.size());
			{
				vector<int> acd(dims.size() + 1);
				acd[0] = 1;
				for (int i = 1; i <= dims.size(); i++) acd[i] = acd[i - 1] * dims[i - i];
				copy(acd.begin() + 1, acd.end(), back_inserter(wAccDims));
			}

			dims = wDims;
			accDims = wAccDims;
		}

		int GetNumDims() const { return dims.size() - 1; }

		int GetDim(int d) {
			assert(d >= 1);
			assert(d <= GetNumDims());
			return dims.at(d);
		}
	};
};

using namespace OneD;

struct Maurer {
	N varN;
	I varI;
	F varF;
};

int main(int argc, char **argv) {
	return EXIT_SUCCESS;
}
