#include <cstdlib>
#include <cassert>
#include <iterator>
#include <vector>
#include <exception>
#include <memory>
#include <iostream>

using namespace std;

class OverflowExc : public ::std::exception {};

/**
= Dimension convention =
dims[1] = x (cols)
dims[2] = y (rows)
== Inside Parse vvLines ==
vvLines    = rows
vvLines[m] = rowM
vvLines[m][n] = rowM, colN
*/

template <typename T>
class VecOb {
public:
	typedef vector<T> vec_t;
	typedef typename vec_t::value_type value_type;
	typedef typename vec_t::reference reference;
	typedef typename vec_t::const_reference const_reference;
	typedef typename vec_t::iterator iterator;
	typedef typename vec_t::const_iterator const_iterator;
	typedef typename vec_t::difference_type difference_type;
	typedef typename vec_t::size_type size_type;
	typedef typename vec_t::pointer pointer;
	vec_t data;

	VecOb() :
		data(1, T()) {}

	VecOb(typename vec_t::size_type maxWanted) :
		data(maxWanted + 1, T()) {}

	VecOb(typename vec_t::size_type maxWanted, const T &a) :
		data(maxWanted + 1, a) {}

	VecOb & operator=(const VecOb &rhs) {
		if (this != &rhs) {
			data = rhs.data;
		}
		return *this;
	}

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

	static VecOb<T> MakeND(int nd, const T &v) {
		VecOb<T> w;

		for (int i = 0; i < nd; i++)
			w.push_back(v);

		return w;
	}

	static VecOb<int> Add(const VecOb<int> &p, const VecOb<int> &r) {
		VecOb<int> w;

		assert(p.size() == r.size());

		for (size_t i = 1; i < p.OnePast(); i++)
			w.push_back(p[i] + r[i]);

		return w;
	}
};

namespace OneD {

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

		VoxelRef(int d) :
			undef(true),
			i(d, 0) {}

	public:

		VoxelRef & operator=(const VoxelRef &rhs) {
			if (this != &rhs) {
				undef = rhs.undef;
				i     = rhs.i;
			}
			return *this;
		}

		static VoxelRef DOf(const VoxelRef &vOther, int dToSet, int newI) {
			VoxelRef w(vOther);

			w.i[dToSet] = newI;

			return w;
		}

		static VoxelRef MakeUndef() {
			return VoxelRef(0);
		}

		static VoxelRef MakeFromVec(const VecOb<int> &v) {
			VoxelRef ret(VoxelRef::MakeUndef());

			ret.undef = false;
			ret.i = v;

			return ret;
		}

		static VoxelRef MakeMergingPrefix(int i1, const VecOb<int> &j) {
			return VoxelRef::MakeFromVec(MergePrefix(i1, j));
		}

		static bool Eq(const VecOb<int> &a, const VecOb<int> &b) {
			assert(a.size() == b.size());

			for (size_t i = 1; i < a.OnePast(); i++)
				if (a[i] != b[i])
					return false;

			return true;
		}

		static bool Eq(const VoxelRef &a, const VoxelRef &b) {
			if (! ((a.undef && b.undef) || (!a.undef && !b.undef)))
				return false;

			return Eq(a.i, b.i);
		}

		static bool Eq(const VoxelRef &a, const VecOb<int> &b) {
			if (a.undef)
				return false;

			return Eq(a.i, b);
		}

		static VoxelRef Add(const VoxelRef &a, const VecOb<int> &b) {
			if (a.undef)
				return a;

			VecOb<int> w;

			assert(a.i.size() == b.size());

			for (size_t i = 1; i < b.OnePast(); i++)
				w.push_back(a.i[i] + b[i]);

			return VoxelRef::MakeFromVec(w);
		}

		static int Euc(const VoxelRef &a) {
			int sum = 0;
			for (size_t i = 1; i < a.i.OnePast(); i++)
				sum += a.i[i] * a.i[i];
			return sum;
		}

		static VoxelRef MinEuc(const VoxelRef &a, const VoxelRef &b) {
			if (a.undef) return b;
			if (b.undef) return a;

			assert(a.i.size() == b.i.size());

			if (Euc(a) < Euc(b)) return a;
			else                 return b;
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
		const int & operator[](const size_t i) const { return dims[i]; }

		size_t TotalAlloc() const {
			int total = 1;

			for (size_t i = 1; i < dims.OnePast(); i++)
				total *= dims[i];

			return total;
		}
	};

	template<typename T>
	struct NStore {
		N dims;
		vector<T> data;

		template<typename DefVal>
		NStore(const VecOb<int> &dms, const DefVal &val) :
			dims(dms),
			data(dims.TotalAlloc(), val) {}

		bool IsInBounds(const VecOb<int> &corVec) const {
			assert(dims.dims.size() == corVec.size());

			for (size_t i = 1; i < corVec.OnePast(); i++)
				if (corVec[i] < 1 || corVec[i] > dims[i])
					return false;

			return true;
		}

		size_t GetIdxOf(const VecOb<int> &corVec) const {
			assert(dims.dims.size() == corVec.size());

			for (auto &i : corVec)
				assert(i > 0);

			size_t pos = 0;
			for (size_t i = 1; i < corVec.OnePast(); i++)
				pos += dims.accDims[i] * (corVec[i] - 1);

			return pos;
		}

		const T & operator[](const VoxelRef &v) const {
			assert(!v.undef);
			return data.at(GetIdxOf(v.i));
		}

		T & operator[](const VoxelRef &v) {
			assert(!v.undef);
			return data.at(GetIdxOf(v.i));
		}

		const T & operator[](const VecOb<int> &v) const {
			return data.at(GetIdxOf(v));
		}

		T & operator[](const VecOb<int> &v) {
			return data.at(GetIdxOf(v));
		}
	};

	struct I : public NStore<int> {
		I(const VecOb<int> &dms) : NStore(dms, 0) {}

		/* d: varying coordinate; r: rest */
		vector<VoxelRef> GetRow(int d, const VoxelRef &r) {
			vector<VoxelRef> vr;

			for (int i = 1; i <= dims[d]; i++)
				vr.push_back(VoxelRef::DOf(r, d, i));

			return vr;
		}
	};

	struct F : public NStore<VoxelRef> {
		F(const VecOb<int> &dms) : NStore(dms, VoxelRef::MakeUndef()) {}
	};

	namespace VecCount {

		void ReinitLower(int d, VecOb<int> *vInOut) {
			assert((int)vInOut->size() >= d - 1);

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

			try {
				Inc_(1, limits, vInOut);
			} catch (const OverflowExc &) {
				isOverflow = true;
			}

			return isOverflow;
		}

		VecOb<int> MakeInitial(const VecOb<int> &limits) {
			return VecOb<int>(limits.size(), 1);
		}

		VecOb<int> MakeLimits(const VecOb<int> dims, int d) {
			/* In the boundary case where d == 1, limits should still produce a sequence,
			such that VecCount::MakeInitial(limit) produces a counter stepping variables i_1...i_(d-1),
			with the result of appending the following: i_1...i_(d-1) :: d :: j_(d+1)...j_k being 'k' of length. */
			VecOb<int> limits;
			copy_n(dims.begin(), d - 1, back_inserter(limits));
			return limits;
		}

	};

	int ChrInt(char c) {
		switch (c) {
		case '0':
			return 0;
		case '1':
			return 1;
		default:
			assert(0);
		}

		return 0xFFFF;
	}
};

using namespace OneD;

namespace Parse {

	typedef VecOb<VecOb<VoxelRef> > vrows_t;
	typedef VecOb<VoxelRef> vrow_t;

	typedef VecOb<VecOb<VecOb<int> > > rows_t;
	typedef VecOb<VecOb<int> > row_t;

	enum CToken {
		TLINE = 0,
		TELT,
		TEND,
	};

	struct D {
		string &s;
		size_t p;

		D(string &s) : s(s), p(0) {}

		bool Eof() {
			return p >= s.size();
		}

		void OptSkipWs() {
			while (!Eof() && (PeekChar() == ' ' || PeekChar() == '\t' || PeekChar() == '\n'))
				p++;
		}

		int PeekChar() {
			assert(!Eof());
			return s[p];
		}

		int GetChar() {
			assert(!Eof());
			return s[p++];
		}

		void ExChar(char c) {
			int w = GetChar();
			assert(w == c);
		}

		CToken CToken() {
			switch (GetChar()) {
			case ';':
				return TLINE;
			case ',':
				return TELT;
			default:
				assert(0);
			};

			return TEND;
		}

		void Unget() {
			assert(p > 0);
			p--;
		}

		int ExNum() {
			assert(!Eof());

			const char *c = s.c_str() + p;
			char *end = NULL;

			long int li = strtol(c, &end, 10);

			if (end == c)
				assert(0);

			p += end - c;

			return li;
		}

		bool TryNum() {
			assert(!Eof());

			const char *c = s.c_str() + p;
			char *end = NULL;

			long int li = strtol(c, &end, 10);

			if (end == c)
				return false;
			else
				return true;
		}

		int ExElt() {
			int w = ExNum();
			OptSkipWs();
			return w;
		}
	};

	void CheckSameSizes(const rows_t &vvLines) {
		size_t *t = nullptr;

		size_t tmp;

		for (auto &i : vvLines)
			for (auto &j : i) {
				if (!t) *(t = &tmp) = j.size();
				assert(*t == j.size());
			}
	}

	void CheckUniformDims(const rows_t &vvLines) {
		size_t numL = vvLines.size();

		for (auto &i : vvLines)
			assert(i.size() == numL);
	}

	void CheckInnerSize(const rows_t &vvLines, int is) {
		for (auto &i : vvLines)
			for (auto &j : i)
				assert(j.size() == is);
	}

	row_t ParseRow(D *d) {
		row_t vLine;

		d->OptSkipWs();

		while (!d->Eof() && d->CToken() == TELT) {
			VecOb<int> iElt;

			while (!d->Eof() && d->TryNum())
				iElt.push_back(d->ExElt());

			vLine.push_back(iElt);
		}

		if (!d->Eof())
			d->Unget();

		return vLine;
	}

	rows_t ParseRows(D *d) {
		rows_t vvLines;

		d->OptSkipWs();

		while (!d->Eof() && d->CToken() == TLINE) {
			vvLines.push_back(ParseRow(d));
		}

		if (!d->Eof())
			d->Unget();

		return vvLines;
	}

	rows_t GetRows(const string &s) {
		string e(s);
		D d(e);

		rows_t vvLines = ParseRows(&d);

		assert(d.Eof());
		CheckSameSizes(vvLines);

		return vvLines;
	}

	VecOb<int> GetRowsDims(const rows_t &vvLines) {
		VecOb<int> w;

		assert(vvLines.size());
		CheckUniformDims(vvLines);

		w.push_back(vvLines[1].size());
		w.push_back(vvLines.size());

		return w;
	}

	VecOb<int> GetUniformDims3(const rows_t &vvLines) {
		assert(vvLines.size());
		const int colLen = vvLines[1].size();
		const int rowsNeeded = colLen * colLen;
		assert(vvLines.size() == rowsNeeded);

		VecOb<int> w;
		w.push_back(colLen);
		w.push_back(colLen);
		w.push_back(colLen);

		return w;
	}
};


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

	static Maurer * Make1DStr(const string &s) {
		VecOb<int> dms(1, s.size());
		Maurer *m = new Maurer(dms);
		for (int i = 1; i <= m->varN[1]; i++)
			m->varI[VoxelRef::MakeFromVec(VecOb<int>(1, i))] = ChrInt(s[i - 1]);
		return m;
	}

	bool Check1D(Parse::row_t row) {
		assert(varN.dims.size() == 1);
		assert(varN[1] == row.size());

		for (int i = 1; i <= varN[1]; i++)
			assert(
			VoxelRef::Eq(
			varF[VecOb<int>(1, i)],
			row[i]));

		return true;
	}

	static Maurer * Make2DStr(const string &s) {
		Parse::rows_t rows = Parse::GetRows(s);
		Parse::CheckInnerSize(rows, 1);

		Maurer *m = new Maurer(Parse::GetRowsDims(rows));

		assert(m->varN.dims.size() == 2);

		for (int i = 1; i <= m->varN[1]; i++)
			for (int j = 1; j <= m->varN[2]; j++) {
				VecOb<int> w; w.push_back(i); w.push_back(j);
				m->varI[w] = rows[i][j][1];
			}

			return m;
	}

	/* FIXME: Not Dist Squared */
	int EucDist(const VoxelRef &a, const VoxelRef &b) const {
		assert(a.i.size() == b.i.size());
		int sum = 0;
		for (size_t i = 1; i < a.i.OnePast(); i++)
			sum += (a.i[i] - b.i[i]) * (a.i[i] - b.i[i]);
		return sum;
	}

	int URdSq(const VoxelRef &u, const VoxelRef &Rd, int d) const {
		assert(u.i.size() == Rd.i.size());
		int sum = 0;
		for (size_t i = 1; i < u.i.OnePast() && i != d; i++)
			sum += (u.i[i] - Rd.i[i]) * (u.i[i] - Rd.i[i]);
		return sum;
	}

	bool RemoveFT(const VoxelRef &u, const VoxelRef &v, const VoxelRef &w, const VoxelRef &Rd, int d) {
		int a = v.i[d] - u.i[d];
		int b = w.i[d] - v.i[d];
		int c = a + b;
		return (c * URdSq(v, Rd, d) - b * URdSq(u, Rd, d) - a * URdSq(w, Rd, d) - a * b * c) > 0;
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
					while (l >= 2 && RemoveFT(gs[l - 1], gs[l], fi, xi, d))
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
				if (varI[w] == 1)
					varF[w] = w;
				else
					varF[w] = u;
			}
		} else {
			for (int id = 1; id <= varN.dims[d]; id++)
				ComputeFT(d - 1, MergePrefix(id, j));
		}

		VecOb<int> limits = VecCount::MakeLimits(varN.dims, d);
		VecOb<int> count  = VecCount::MakeInitial(limits);

		do {
			VoronoiFT(d, count, j);
		} while (!VecCount::Inc(limits, &count));
	}

	void Start() {
		ComputeFT(varN.dims.size(), VecOb<int>(0));
	}
};

namespace B84d {

	typedef NStore<VoxelRef> F;

	/* FIXME: (x,y) or (y,x) ? Borgefors84::Algorithm 1. uses (y,x)
	Decided on use of (x,y) for consistency */
	struct Mask {
		size_t n;

		VecOb<VecOb<int> > pos;
		VecOb<VecOb<int> > add;

		/* Multieval */
#define MASK_MAKE(n, p, a) (Make((n), (p), sizeof((p))/sizeof((*p)), (a), sizeof((a))/sizeof((*a))))

		static Mask Make(size_t n, int p[], size_t pLen, int a[], size_t aLen) {
			Mask r;

			r.n = n;

			assert(pLen == aLen);
			assert(pLen % n == 0);

			for (size_t i = 0; i < pLen; i += n) {
				VecOb<int> nPos;
				VecOb<int> nAdd;

				for (size_t j = 0; j < n; j++)
					nPos.push_back(p[i+j]);
				for (size_t j = 0; j < n; j++)
					nAdd.push_back(a[i+j]);

				r.pos.push_back(nPos);
				r.add.push_back(nAdd);
			}

			return r;
		}

		static Mask Make2D_FF() {
			int p[] = {
				-1, -1, 0, -1, 1, -1,
				-1,  0, 0,  0,
			};

			int a[] = {
				1, 1, 0, 1, 1, 1,
				1, 0, 0, 0,
			};

			return MASK_MAKE(2, p, a);
		}

		static Mask Make2D_FB() {
			int p[] = {
				0, 0, 1, 0,
			};

			int a[] = {
				0, 0, 1, 0,
			};

			return MASK_MAKE(2, p, a);
		}

		static Mask Make2D_BF() {
			int p[] = {
				0, 0, 1, 0,
				-1, 1, 0, 1, 1, 1,
			};

			int a[] = {
				0, 0, 1, 0,
				1, 1, 0, 1, 1, 1,
			};

			return MASK_MAKE(2, p, a);
		}

		static Mask Make2D_BB() {
			int p[] = {
				-1, 0, 0, 0,
			};

			int a[] = {
				1, 0, 0, 0,
			};

			return MASK_MAKE(2, p, a);
		}

		size_t size() const {
			return pos.size();
		}

		VoxelRef LowestFAround(const F &varF, const VecOb<int> &p) const {
			VecOb<VoxelRef> vox(0, VoxelRef::MakeUndef());

			for (size_t i = 1; i <= size(); i++) {
				VecOb<int> candidate = VecOb<int>::Add(p, pos[i]);
				if (varF.IsInBounds(candidate))
					vox.push_back(VoxelRef::Add(varF[candidate], add[i]));
			}

			assert(vox.size());

			VoxelRef cur = vox[1];

			for (auto &i : vox)
				cur = VoxelRef::MinEuc(cur, i);

			return cur;
		}
	};

	struct DEuc {
		N varN;
		I varI;
		F varF;

	private:
		DEuc(const VecOb<int> &dms) :
			varN(dms),
			varI(dms),
			varF(dms, VoxelRef::MakeFromVec(VecOb<int>::MakeND(dms.size(), 0))) {}

	public:

		void InitIFElt(const VecOb<int> &w, int rowsVal) {
			varI[w] = rowsVal;

			/* Borgefors84::2:'zero for feature elements and infinity otherwise' */
			switch (rowsVal) {
			case 0:
				varF[w] = VoxelRef::MakeUndef();
				break;
			case 1:
				varF[w] = VoxelRef::MakeFromVec(VecOb<int>::MakeND(varN.dims.size(), 0));
				break;
			default:
				assert(0);
			}
		}

		static DEuc * Make2DStr(const string &s) {
			Parse::rows_t rows = Parse::GetRows(s);
			Parse::CheckInnerSize(rows, 1);

			DEuc *m = new DEuc(Parse::GetRowsDims(rows));

			assert(m->varN.dims.size() == 2);

			for (int r = 1; r <= m->varN[2]; r++)
				for (int c = 1; c <= m->varN[1]; c++) {
					VecOb<int> w; w.push_back(c); w.push_back(r);

					m->InitIFElt(w, rows[r][c][1]);
				}

				return m;
		}

		static DEuc * Make3DUniformStr(const string &s) {
			Parse::rows_t rows = Parse::GetRows(s);
			Parse::CheckInnerSize(rows, 1);

			DEuc *m = new DEuc(Parse::GetUniformDims3(rows));

			for (int s = 1; s <= m->varN[3]; s++)
				for (int r = 1; r <= m->varN[2]; r++)
					for (int c = 1; c <= m->varN[1]; c++) {
						VecOb<int> w; w.push_back(c); w.push_back(r); w.push_back(s);

						m->InitIFElt(w, rows[(m->varN[2] * (s-1)) + r][c][1]);
					}

					return m;
		}

		void Start() {
		}

		void Start2D() {
			assert(varN.dims.size() == 2);

			Mask maskFF = Mask::Make2D_FF();
			Mask maskFB = Mask::Make2D_FB();

			Mask maskBF = Mask::Make2D_BF();
			Mask maskBB = Mask::Make2D_BB();

			/* Forward pass */
			for (int r = 2; r <= varN.dims[2]; r++) {
				/* L->R */
				for (int c = 2; c <= varN.dims[1]; c++) {
					VecOb<int> w; w.push_back(c); w.push_back(r);
					varF[w] = maskFF.LowestFAround(varF, w);
				}

				/* R->L */
				for (int c = varN.dims[1] - 1; c >= 1; c--) {
					VecOb<int> w; w.push_back(c); w.push_back(r);
					varF[w] = maskFB.LowestFAround(varF, w);
				}
			}

			/* Backward pass */
			for (int r = varN.dims[2]; r >= 1; r--) {
				/* R->L */
				for (int c = varN.dims[1]; c >= 1; c--) {
					VecOb<int> w; w.push_back(c); w.push_back(r);
					varF[w] = maskBF.LowestFAround(varF, w);
				}

				/* L->R */
				for (int c = 1; c <= varN.dims[1]; c++) {
					VecOb<int> w; w.push_back(c); w.push_back(r);
					varF[w] = maskBB.LowestFAround(varF, w);
				}
			}
		}

		void Start3D() {

		}
	};

};

template<typename T>
void PFELT2(const T &varF, const VecOb<int> &a) {
	if (varF[a].undef) cout << " ," << "X X";
	else               cout << " ," << varF[a].i[1] << " " << varF[a].i[2];
}

template<typename T>
void PF2(const T &varF) {
	assert(varF.dims.dims.size() == 2);
	for (int r = 1; r <= varF.dims[2]; r++) {
		cout << ";";
		for (int c = 1; c <= varF.dims[1]; c++) {
			VecOb<int> w; w.push_back(c); w.push_back(r);
			PFELT2(varF, w);
		}
		cout << endl;
	}
	cout << endl;
}

void PM2(const Maurer &m) {
	PF2(m.varF);
}

void PD2(const B84d::DEuc &m) {
	PF2(m.varF);
}

void T1DStr() {
	shared_ptr<Maurer> m(Maurer::Make1DStr("10010"));
	m->Start();
	Parse::rows_t rows = Parse::GetRows(";,1,1,4,4,4");
	assert(rows.size() == 1);
	m->Check1D(rows[1]);
}

void T2DStr() {
	Maurer *pm;
	shared_ptr<Maurer> m((pm = Maurer::Make2DStr("; ,1 ,0 ,0 ; ,0 ,1 ,0 ; ,0 ,0 ,0")));
	m->Start();
}

void TE2DStr() {
	using namespace B84d;
	DEuc *pm;
	shared_ptr<DEuc> m((pm = DEuc::Make2DStr("; ,1 ,0 ,0 ; ,0 ,1 ,0 ; ,0 ,0 ,0")));
	m->Start2D();
}

void TE3DStr() {
	using namespace B84d;
	DEuc *pm;
	shared_ptr<DEuc> m((pm = DEuc::Make3DUniformStr(
		"; ,1 ,0 ,0 ; ,0 ,1 ,0 ; ,0 ,0 ,0"
		"; ,0 ,0 ,0 ; ,0 ,0 ,0 ; ,0 ,0 ,0"
		"; ,0 ,0 ,0 ; ,0 ,0 ,0 ; ,0 ,0 ,0")));
	m->Start3D();
}

int main(int argc, char **argv) {
	//T1DStr();
	T2DStr();
	TE2DStr();
	TE3DStr();
	return EXIT_SUCCESS;
}
