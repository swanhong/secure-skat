package ring

import (
	"fmt"
	"unsafe"
)

type NumberTheoreticTransformer3N struct {
	numberTheoreticTransformerBase

	N    int    // Transform size
	p    uint64 // Prime modulus
	a, b int    // N = 2^a * 3^b factorization

	zetasMont       []uint64
	omegaMont       uint64
	NInvMont        uint64
	zMinusZ5InvMont uint64
}

// NewNumberTheoreticTransformer3N creates a 3N-NTT transformer.
func NewNumberTheoreticTransformer3N(r *SubRing, _ int) NumberTheoreticTransformer {
	q := r.Modulus
	n := r.N
	bred := r.BRedConstant
	mred := r.MRedConstant

	a, b := r.a, r.b
	if n <= 0 || a <= 0 {
		panic(fmt.Sprintf("3N-NTT requires initialized SubRing factors a,b for N=%d", n))
	}

	var w uint64
	if r.NTTTable != nil {
		w = r.PrimitiveRoot
	}
	if w == 0 || OrderMod(w, q) != uint64(3*n) {
		var err error
		w, err = FindPrimitiveRootOfUnity(q, uint64(3*n))
		if err != nil {
			panic(fmt.Sprintf("failed to find primitive 3N-th root: %v", err))
		}
	}

	omega := ModExp(w, uint64(n), q)
	z := ModExp(w, uint64(n/2), q)

	NInvMont := MForm(ModExp(uint64(n), q-2, q), q, bred)
	z5 := ModExp(z, 5, q) // z^5 mod q
	zMinusZ5 := z + q - z5
	zetasMont := genCyclotomicTwiddleFactors(n, a, b, w, q, bred)
	omegaMont := MForm(omega, q, bred)
	zMinusZ5InvMont := MForm(ModExp(zMinusZ5, q-2, q), q, bred)

	nttTable := r.NTTTable
	if nttTable == nil {
		nttTable = &NTTTable{
			NthRoot:       uint64(3 * n),
			PrimitiveRoot: w,
			RootsForward:  make([]uint64, n),
			RootsBackward: make([]uint64, n),
			NInv:          NInvMont,
		}
	}

	ntt := &NumberTheoreticTransformer3N{
		numberTheoreticTransformerBase: numberTheoreticTransformerBase{
			N:            n,
			Modulus:      q,
			MRedConstant: mred,
			BRedConstant: bred,
			NTTTable:     nttTable,
		},
		N:               n,
		p:               q,
		a:               a,
		b:               b,
		zetasMont:       zetasMont,
		omegaMont:       omegaMont,
		NInvMont:        NInvMont,
		zMinusZ5InvMont: zMinusZ5InvMont,
	}

	return ntt
}

func genCyclotomicTwiddleFactors(N, radix2, radix3 int, root, q uint64, bred [2]uint64) []uint64 {
	level := radix2 + radix3
	tree := make([][]uint64, level+1)
	for i := range tree {
		tree[i] = make([]uint64, N)
	}

	order := uint64(3 * N)
	tree[0][0] = order
	tree[1][0] = order / 6
	tree[1][1] = 5 * order / 6

	zetasMont := []uint64{
		MForm(1, q, bred),
		MForm(ModExp(root, tree[1][0], q), q, bred),
	}
	width := 2
	for l := 1; l <= radix3; l++ {
		for i := 0; i < width; i++ {
			base := tree[l][i] / 3
			tree[l+1][3*i] = base
			tree[l+1][3*i+1] = base + order/3
			tree[l+1][3*i+2] = base + 2*order/3
			zetasMont = append(zetasMont, MForm(ModExp(root, base, q), q, bred))
			zetasMont = append(zetasMont, MForm(ModExp(root, 2*base, q), q, bred))
		}
		width *= 3
	}

	for l := radix3 + 1; l < level; l++ {
		for i := 0; i < width; i++ {
			base := tree[l][i] / 2
			tree[l+1][2*i] = base
			tree[l+1][2*i+1] = base + order/2
			zetasMont = append(zetasMont, MForm(ModExp(root, base, q), q, bred))
		}
		width *= 2
	}

	return zetasMont
}

func butterfly3(U, V, W, Psi1, Psi2, twoQ, fourQ, Q, omegaMont, MRedConstant uint64) (uint64, uint64, uint64) {
	if U >= fourQ {
		U -= fourQ
	}

	T1 := MRedLazy(V, Psi1, Q, MRedConstant) // z1 * p2[i+step]
	T2 := MRedLazy(W, Psi2, Q, MRedConstant) // z2 * p2[i+2*step]

	T3 := T1 + twoQ - T2                          // t1 - t2
	T3 = MRedLazy(T3, omegaMont, Q, MRedConstant) // ω*(t1 - t2)

	W = U + fourQ - T1 - T3 // b - t1 - t3
	V = U + twoQ - T2 + T3  // b - t2 + t3
	U = U + T1 + T2         // b + t1 + t2

	return U, V, W
}

func invbutterfly3(U, V, W, Psi1, Psi2, twoQ, fourQ, Q, omegaMont, MRedConstant uint64) (X, Y, Z uint64) {

	X = addLazy(addLazy(U, V, twoQ), W, twoQ)

	T := V + twoQ - U
	T = MRedLazy(T, omegaMont, Q, MRedConstant)

	Y = W + twoQ - U + T
	Y = MRedLazy(Y, Psi1, Q, MRedConstant)

	Z = W + fourQ - V - T
	Z = MRedLazy(Z, Psi2, Q, MRedConstant)

	return
}

func (ntt *NumberTheoreticTransformer3N) Forward(p1, p2 []uint64) {
	n := ntt.N
	if len(p1) < n || len(p2) < n {
		panic(fmt.Sprintf("Forward: len(p1)=%d len(p2)=%d < N=%d", len(p1), len(p2), n))
	}

	ntt.FactorizedDFT(p1, p2)
	reducevec(p2, p2, ntt.p, ntt.BRedConstant)
}

// FactorizedDFT computes the fast O(N log N) cyclotomic polynomial DFT
// (Montgomery-aware: roots in Montgomery form, coeffs standard; all mult via MRed/MRedLazy)
func (ntt *NumberTheoreticTransformer3N) FactorizedDFT(p1, p2 []uint64) {

	// Sanity check
	if len(p2) < MinimumRingDegreeForLoopUnrolledNTT {
		panic(fmt.Sprintf("unsafe call of FactorizedDFT: receiver len(p2)=%d < %d", len(p2), MinimumRingDegreeForLoopUnrolledNTT))
	}

	var t int
	var F, V uint64

	fourQ := ntt.p << 2
	twoQ := ntt.p << 1

	// 1) Special initial cyclotomic layer
	t = ntt.N >> 1
	F = ntt.zetasMont[1]

	for jx, jy := 0, t; jx < t; jx, jy = jx+8, jy+8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
		xin := (*[8]uint64)(unsafe.Pointer(&p1[jx]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
		yin := (*[8]uint64)(unsafe.Pointer(&p1[jy]))

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
		xout := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
		yout := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

		V = MRedLazy(yin[0], F, ntt.p, ntt.MRedConstant)
		xout[0], yout[0] = xin[0]+V, xin[0]+yin[0]+twoQ-V

		V = MRedLazy(yin[1], F, ntt.p, ntt.MRedConstant)
		xout[1], yout[1] = xin[1]+V, xin[1]+yin[1]+twoQ-V

		V = MRedLazy(yin[2], F, ntt.p, ntt.MRedConstant)
		xout[2], yout[2] = xin[2]+V, xin[2]+yin[2]+twoQ-V

		V = MRedLazy(yin[3], F, ntt.p, ntt.MRedConstant)
		xout[3], yout[3] = xin[3]+V, xin[3]+yin[3]+twoQ-V

		V = MRedLazy(yin[4], F, ntt.p, ntt.MRedConstant)
		xout[4], yout[4] = xin[4]+V, xin[4]+yin[4]+twoQ-V

		V = MRedLazy(yin[5], F, ntt.p, ntt.MRedConstant)
		xout[5], yout[5] = xin[5]+V, xin[5]+yin[5]+twoQ-V

		V = MRedLazy(yin[6], F, ntt.p, ntt.MRedConstant)
		xout[6], yout[6] = xin[6]+V, xin[6]+yin[6]+twoQ-V

		V = MRedLazy(yin[7], F, ntt.p, ntt.MRedConstant)
		xout[7], yout[7] = xin[7]+V, xin[7]+yin[7]+twoQ-V
	}

	zetaIdx := 2

	// 2) Radix-3 layers
	t = ntt.N / 6
	radix3 := ntt.b

	for radix3 >= 1 {
		for start := 0; start < ntt.N; start += 3 * t {
			z1 := ntt.zetasMont[zetaIdx]
			z2 := ntt.zetasMont[zetaIdx+1]
			zetaIdx += 2

			for jx, jy, jz := start, start+t, start+2*t; jx < start+t; jx, jy, jz = jx+8, jy+8, jz+8 {

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
				x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
				y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
				z := (*[8]uint64)(unsafe.Pointer(&p2[jz]))

				x[0], y[0], z[0] = butterfly3(x[0], y[0], z[0], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[1], y[1], z[1] = butterfly3(x[1], y[1], z[1], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[2], y[2], z[2] = butterfly3(x[2], y[2], z[2], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[3], y[3], z[3] = butterfly3(x[3], y[3], z[3], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[4], y[4], z[4] = butterfly3(x[4], y[4], z[4], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[5], y[5], z[5] = butterfly3(x[5], y[5], z[5], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[6], y[6], z[6] = butterfly3(x[6], y[6], z[6], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[7], y[7], z[7] = butterfly3(x[7], y[7], z[7], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
			}
		}
		radix3 = radix3 - 1
		t /= 3
	}

	var reduce bool

	// 3) Radix-2 layers
	t = 1 << (ntt.a - 2)

	for m := 1; m < ntt.a; m++ {

		reduce = m&1 != 0

		if t >= 8 {
			for start := 0; start < ntt.N; start += t << 1 {
				F = ntt.zetasMont[zetaIdx]
				zetaIdx++

				if reduce {

					for jx, jy := start, start+t; jx < start+t; jx, jy = jx+8, jy+8 {

						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						x[0], y[0] = butterfly(x[0], y[0], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
						x[1], y[1] = butterfly(x[1], y[1], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
						x[2], y[2] = butterfly(x[2], y[2], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
						x[3], y[3] = butterfly(x[3], y[3], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
						x[4], y[4] = butterfly(x[4], y[4], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
						x[5], y[5] = butterfly(x[5], y[5], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
						x[6], y[6] = butterfly(x[6], y[6], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
						x[7], y[7] = butterfly(x[7], y[7], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					}
				} else {

					for jx, jy := start, start+t; jx < start+t; jx, jy = jx+8, jy+8 {

						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						V = MRedLazy(y[0], F, ntt.p, ntt.MRedConstant)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedLazy(y[1], F, ntt.p, ntt.MRedConstant)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedLazy(y[2], F, ntt.p, ntt.MRedConstant)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedLazy(y[3], F, ntt.p, ntt.MRedConstant)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedLazy(y[4], F, ntt.p, ntt.MRedConstant)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedLazy(y[5], F, ntt.p, ntt.MRedConstant)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedLazy(y[6], F, ntt.p, ntt.MRedConstant)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedLazy(y[7], F, ntt.p, ntt.MRedConstant)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}

			}

		} else if t >= 4 {
			if reduce {

				for jx := 0; jx < ntt.N; jx, zetaIdx = jx+16, zetaIdx+2 {

					z := (*[2]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[jx]))

					x[0], x[4] = butterfly(x[0], x[4], z[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[1], x[5] = butterfly(x[1], x[5], z[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[2], x[6] = butterfly(x[2], x[6], z[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[3], x[7] = butterfly(x[3], x[7], z[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[8], x[12] = butterfly(x[8], x[12], z[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[9], x[13] = butterfly(x[9], x[13], z[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[10], x[14] = butterfly(x[10], x[14], z[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[11], x[15] = butterfly(x[11], x[15], z[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				}

				reduce = false

			} else {

				for jx := 0; jx < ntt.N; jx, zetaIdx = jx+16, zetaIdx+2 {

					z := (*[2]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[jx]))

					V = MRedLazy(x[4], z[0], ntt.p, ntt.MRedConstant)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[5], z[0], ntt.p, ntt.MRedConstant)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], z[0], ntt.p, ntt.MRedConstant)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[7], z[0], ntt.p, ntt.MRedConstant)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedLazy(x[12], z[1], ntt.p, ntt.MRedConstant)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[13], z[1], ntt.p, ntt.MRedConstant)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], z[1], ntt.p, ntt.MRedConstant)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[15], z[1], ntt.p, ntt.MRedConstant)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V
				}
			}

		} else if t >= 2 {
			if reduce {
				for jx := 0; jx < ntt.N; jx, zetaIdx = jx+16, zetaIdx+4 {

					z := (*[4]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[jx]))

					x[0], x[2] = butterfly(x[0], x[2], z[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[1], x[3] = butterfly(x[1], x[3], z[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[4], x[6] = butterfly(x[4], x[6], z[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[5], x[7] = butterfly(x[5], x[7], z[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[8], x[10] = butterfly(x[8], x[10], z[2], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[9], x[11] = butterfly(x[9], x[11], z[2], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[12], x[14] = butterfly(x[12], x[14], z[3], twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[13], x[15] = butterfly(x[13], x[15], z[3], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				}

				reduce = false

			} else {
				for jx := 0; jx < ntt.N; jx, zetaIdx = jx+16, zetaIdx+4 {

					z := (*[4]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[jx]))

					V = MRedLazy(x[2], z[0], ntt.p, ntt.MRedConstant)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], z[0], ntt.p, ntt.MRedConstant)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], z[1], ntt.p, ntt.MRedConstant)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], z[1], ntt.p, ntt.MRedConstant)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedLazy(x[10], z[2], ntt.p, ntt.MRedConstant)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], z[2], ntt.p, ntt.MRedConstant)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], z[3], ntt.p, ntt.MRedConstant)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], z[3], ntt.p, ntt.MRedConstant)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}
		} else {
			for jx := 0; jx < ntt.N; jx, zetaIdx = jx+16, zetaIdx+8 {

				z := (*[8]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%16 != 0 */
				x := (*[16]uint64)(unsafe.Pointer(&p2[jx]))

				x[0], x[1] = butterfly(x[0], x[1], z[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				x[2], x[3] = butterfly(x[2], x[3], z[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				x[4], x[5] = butterfly(x[4], x[5], z[2], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				x[6], x[7] = butterfly(x[6], x[7], z[3], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				x[8], x[9] = butterfly(x[8], x[9], z[4], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				x[10], x[11] = butterfly(x[10], x[11], z[5], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				x[12], x[13] = butterfly(x[12], x[13], z[6], twoQ, fourQ, ntt.p, ntt.MRedConstant)
				x[14], x[15] = butterfly(x[14], x[15], z[7], twoQ, fourQ, ntt.p, ntt.MRedConstant)
			}
		}
		t >>= 1
	}
}

// FactorizedIDFT computes the fast O(N log N) inverse cyclotomic DFT
// (Montgomery-aware; final 1/N, 2/N scaling done here)
func (ntt *NumberTheoreticTransformer3N) FactorizedIDFT(p1, p2 []uint64) {

	// Sanity check
	if len(p2) < MinimumRingDegreeForLoopUnrolledNTT {
		panic(fmt.Sprintf("unsafe call of FactorizedIDFT: receiver len(p2)=%d < %d", len(p2), MinimumRingDegreeForLoopUnrolledNTT))
	}

	var T1, T2 uint64
	var F uint64

	twoQ := ntt.p << 1
	fourQ := ntt.p << 2

	// 1) Inverse radix-2 layers
	step := 1
	maxStep := 1 << (ntt.a - 2)

	zetaIdx := len(ntt.zetasMont) - 8

	for j := 0; j < ntt.N; j, zetaIdx = j+16, zetaIdx-8 {
		psi := (*[8]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))
		xin := (*[16]uint64)(unsafe.Pointer(&p1[j]))
		xout := (*[16]uint64)(unsafe.Pointer(&p2[j]))

		xout[0], xout[1] = invbutterfly(xin[1], xin[0], psi[7], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		xout[2], xout[3] = invbutterfly(xin[3], xin[2], psi[6], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		xout[4], xout[5] = invbutterfly(xin[5], xin[4], psi[5], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		xout[6], xout[7] = invbutterfly(xin[7], xin[6], psi[4], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		xout[8], xout[9] = invbutterfly(xin[9], xin[8], psi[3], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		xout[10], xout[11] = invbutterfly(xin[11], xin[10], psi[2], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		xout[12], xout[13] = invbutterfly(xin[13], xin[12], psi[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		xout[14], xout[15] = invbutterfly(xin[15], xin[14], psi[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
	}

	step <<= 1

	zetaIdx = zetaIdx + 4

	for j := 0; j < ntt.N; j, zetaIdx = j+16, zetaIdx-4 {
		psi := (*[4]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))
		x := (*[16]uint64)(unsafe.Pointer(&p2[j]))

		x[0], x[2] = invbutterfly(x[2], x[0], psi[3], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[1], x[3] = invbutterfly(x[3], x[1], psi[3], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[4], x[6] = invbutterfly(x[6], x[4], psi[2], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[5], x[7] = invbutterfly(x[7], x[5], psi[2], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[8], x[10] = invbutterfly(x[10], x[8], psi[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[9], x[11] = invbutterfly(x[11], x[9], psi[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[12], x[14] = invbutterfly(x[14], x[12], psi[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[13], x[15] = invbutterfly(x[15], x[13], psi[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
	}

	step <<= 1

	zetaIdx = zetaIdx + 2

	for j := 0; j < ntt.N; j, zetaIdx = j+16, zetaIdx-2 {
		psi := (*[2]uint64)(unsafe.Pointer(&ntt.zetasMont[zetaIdx]))
		x := (*[16]uint64)(unsafe.Pointer(&p2[j]))

		x[0], x[4] = invbutterfly(x[4], x[0], psi[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[1], x[5] = invbutterfly(x[5], x[1], psi[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[2], x[6] = invbutterfly(x[6], x[2], psi[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[3], x[7] = invbutterfly(x[7], x[3], psi[1], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[8], x[12] = invbutterfly(x[12], x[8], psi[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[9], x[13] = invbutterfly(x[13], x[9], psi[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[10], x[14] = invbutterfly(x[14], x[10], psi[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
		x[11], x[15] = invbutterfly(x[15], x[11], psi[0], twoQ, fourQ, ntt.p, ntt.MRedConstant)
	}

	step <<= 1

	zetaIdx = zetaIdx + 1

	for step <= maxStep {

		if step >= 8 {
			for start := 0; start < ntt.N; start, zetaIdx = start+2*step, zetaIdx-1 {
				F = ntt.zetasMont[zetaIdx]

				for jx, jy := start, start+step; jx < start+step; jx, jy = jx+8, jy+8 {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
					x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
					y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

					x[0], y[0] = invbutterfly(y[0], x[0], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[1], y[1] = invbutterfly(y[1], x[1], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[2], y[2] = invbutterfly(y[2], x[2], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[3], y[3] = invbutterfly(y[3], x[3], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[4], y[4] = invbutterfly(y[4], x[4], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[5], y[5] = invbutterfly(y[5], x[5], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[6], y[6] = invbutterfly(y[6], x[6], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
					x[7], y[7] = invbutterfly(y[7], x[7], F, twoQ, fourQ, ntt.p, ntt.MRedConstant)
				}
			}
		}

		step <<= 1
	}

	// 2) Inverse radix-3 layers
	for step <= ntt.N/6 {
		for start := 0; start < ntt.N; start += 3 * step {
			z2 := ntt.zetasMont[zetaIdx]
			zetaIdx--
			z1 := ntt.zetasMont[zetaIdx]
			zetaIdx--

			for i := start; i < start+step; i = i + 8 {
				x := (*[8]uint64)(unsafe.Pointer(&p2[i]))
				y := (*[8]uint64)(unsafe.Pointer(&p2[i+step]))
				z := (*[8]uint64)(unsafe.Pointer(&p2[i+2*step]))

				x[0], y[0], z[0] = invbutterfly3(x[0], y[0], z[0], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[1], y[1], z[1] = invbutterfly3(x[1], y[1], z[1], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[2], y[2], z[2] = invbutterfly3(x[2], y[2], z[2], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[3], y[3], z[3] = invbutterfly3(x[3], y[3], z[3], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[4], y[4], z[4] = invbutterfly3(x[4], y[4], z[4], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[5], y[5], z[5] = invbutterfly3(x[5], y[5], z[5], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[6], y[6], z[6] = invbutterfly3(x[6], y[6], z[6], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
				x[7], y[7], z[7] = invbutterfly3(x[7], y[7], z[7], z1, z2, twoQ, fourQ, ntt.p, ntt.omegaMont, ntt.MRedConstant)
			}
		}
		step *= 3
	}

	// 3) Final cyclotomic layer (inverse), with (1/N) and (2/N) scaling
	half := ntt.N >> 1
	for i := 0; i < half; i = i + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p2[i]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[i+half]))

		T1 = x[0] + y[0]
		T2 = x[0] + twoQ - y[0]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[0] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[0] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T1 = x[1] + y[1]
		T2 = x[1] + twoQ - y[1]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[1] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[1] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T1 = x[2] + y[2]
		T2 = x[2] + twoQ - y[2]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[2] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[2] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T1 = x[3] + y[3]
		T2 = x[3] + twoQ - y[3]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[3] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[3] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T1 = x[4] + y[4]
		T2 = x[4] + twoQ - y[4]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[4] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[4] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T1 = x[5] + y[5]
		T2 = x[5] + twoQ - y[5]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[5] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[5] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T1 = x[6] + y[6]
		T2 = x[6] + twoQ - y[6]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[6] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[6] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T1 = x[7] + y[7]
		T2 = x[7] + twoQ - y[7]
		T2 = MRedLazy(T2, ntt.zMinusZ5InvMont, ntt.p, ntt.MRedConstant)

		T1 = T1 + twoQ - T2
		x[7] = MRed(T1, ntt.NInvMont, ntt.p, ntt.MRedConstant)

		T2 = T2 + T2
		y[7] = MRed(T2, ntt.NInvMont, ntt.p, ntt.MRedConstant)
	}
}

func (ntt *NumberTheoreticTransformer3N) Backward(p1, p2 []uint64) {
	n := ntt.N
	if len(p1) < n || len(p2) < n {
		panic(fmt.Sprintf("Backward: len(p1)=%d len(p2)=%d < N=%d", len(p1), len(p2), n))
	}

	ntt.FactorizedIDFT(p1, p2)
	reducevec(p2, p2, ntt.p, ntt.BRedConstant)
}

// ForwardLazy is a lazy version of Forward that doesn't perform modular reduction
func (ntt *NumberTheoreticTransformer3N) ForwardLazy(p1, p2 []uint64) {
	ntt.Forward(p1, p2)
}

// BackwardLazy is a lazy version of Backward that doesn't perform modular reduction
func (ntt *NumberTheoreticTransformer3N) BackwardLazy(p1, p2 []uint64) {
	ntt.Backward(p1, p2)
}

// ---- lazy helpers (keep values in [0, 2Q)) ----
func addLazy(a, b, twoQ uint64) uint64 {
	s := a + b
	if s >= twoQ {
		s -= twoQ
	}
	return s
}
