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
	struct IData {
		vector<int> dims;
		vector<int> data;

		vector<int> accDims;

		IData(const vector<int> dms)
		{
			assert(dims.size());

			vector<int> wDims = dms;

			vector<int> wData;
			{
				int total = 1;
				for (int i = 0; i < dims.size(); i++) total *= dims[i];
				wData = vector<int>(total, 0);
			}

			vector<int> wAccDims(dims.size());
			{
				vector<int> acd(dims.size() + 1);
				acd[0] = 1;
				for (int i = 1; i <= dims.size(); i++) acd[i] = acd[i - 1] * dims[i - i];
				copy(acd.begin() + 1, acd.end(), back_inserter(wAccDims));
			}

			dims    = wDims;
			data    = wData;
			accDims = wAccDims;
		}

		void Check() {
			assert(dims.size() >= 1);
			assert(dims.size() == accDims.size());
		}

		int GetNumDims() const { return dims.size(); }

		int GetDim(int d) {
			assert(d >= 1);
			assert(d <= GetNumDims());
			return dims.at(d - 1);
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
	};

	struct I {
		IData data;

		I(const vector<int> dms) : data(dms) {}

		int GetNumDims() const { return data.GetNumDims(); }

		int GetDim(int d) { return data.GetDim(d); }

		/* d: varying coordinate; r: rest */
		vector<VoxelRef> GetRow(int d, const VoxelRef &r) {
			vector<VoxelRef> vr;

			for (int i = 1; i <= GetDim(d); i++)
				vr.push_back(VoxelRef::DOf(r, d, i));

			return vr;
		}
	};

	struct SubImage {
		SubImage(const I &img) {}
	};
};

int main(int argc, char **argv) {
	return EXIT_SUCCESS;
}
