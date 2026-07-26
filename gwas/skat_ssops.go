package gwas

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// Fixed-point secret-sharing and shape helpers for secure SKAT.

// ssMul returns TruncVec(a⊙b) at fracBits — secret×secret elementwise multiply.
func (ast *AssocTest) ssMul(a, b mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	return mpcObj.TruncVec(mpcObj.SSMultElemVec(a, b), mpcObj.GetDataBits(), mpcObj.GetFracBits())
}

// ssSquare returns TruncVec(a⊙a) at fracBits.
func (ast *AssocTest) ssSquare(a mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	return mpcObj.TruncVec(mpcObj.SSSquareElemVec(a), mpcObj.GetDataBits(), mpcObj.GetFracBits())
}

// ssPMul returns TruncVec(a·cf) — secret vector times a public constant.
func (ast *AssocTest) ssPMul(a mpc_core.RVec, cf float64) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	cfE := mpcObj.GetRType().FromFloat64(cf, fb)
	out := mpc_core.RMultConstVec(cfE, a)
	return mpcObj.TruncVec(out, db, fb)
}

// ssDot returns Σᵢ TruncVec(a⊙b)ᵢ — secret dot product reduced to one scalar.
func (ast *AssocTest) ssDot(a, b mpc_core.RVec) mpc_core.RElem {
	mpcObj := ast.general.mpcObj[0]
	p := mpcObj.TruncVec(mpcObj.SSMultElemVec(a, b), mpcObj.GetDataBits(), mpcObj.GetFracBits())
	acc := mpcObj.GetRType().Zero()
	for i := range p {
		acc = acc.Add(p[i])
	}
	return acc
}

// ssMatVec returns the untruncated secret-shared product M·v.
func (ast *AssocTest) ssMatVec(M mpc_core.RMat, v mpc_core.RVec) mpc_core.RVec {
	return col0(ast.general.mpcObj[0].SSMultMat(M, asCol(v)))
}

// hubVec returns an additive share of the public constant cnst (hub holds it at fracBits, others 0).
func (ast *AssocTest) hubVec(cnst float64, n int) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	e := mpcObj.GetRType().Zero()
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		e = mpcObj.GetRType().FromFloat64(cnst, mpcObj.GetFracBits())
	}
	return mpc_core.InitRVec(e, n)
}

// asCol wraps v as an n×1 column RMat; col0 extracts column 0. Pure reshape, no arithmetic/truncation.
func asCol(v mpc_core.RVec) mpc_core.RMat {
	out := make(mpc_core.RMat, len(v))
	for i := range v {
		out[i] = mpc_core.RVec{v[i]}
	}
	return out
}

func col0(M mpc_core.RMat) mpc_core.RVec {
	out := make(mpc_core.RVec, len(M))
	for i := range M {
		out[i] = M[i][0]
	}
	return out
}

// firstCt returns element 0 of a CipherVector, or nil when empty (the pid-0 aggregate case).
func firstCt(cv crypto.CipherVector) *rlwe.Ciphertext {
	if len(cv) > 0 {
		return cv[0]
	}
	return nil
}
