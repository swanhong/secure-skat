package gwas

import mpc_core "github.com/hhcho/mpc-core"

func matrixDims(a mpc_core.RMat) (rows, columns int) {
	rows = len(a)
	if rows == 0 {
		return
	}
	columns = len(a[0])
	for i := 1; i < rows; i++ {
		if len(a[i]) != columns {
			panic("ragged secret-shared matrix")
		}
	}
	return
}

func transposeShares(rtype mpc_core.RElem, a mpc_core.RMat) mpc_core.RMat {
	rows, columns := matrixDims(a)
	out := mpc_core.InitRMat(rtype.Zero(), columns, rows)
	for i := 0; i < rows; i++ {
		for j := 0; j < columns; j++ {
			out[j][i] = a[i][j]
		}
	}
	return out
}

func matrixFromFlat(values mpc_core.RVec, rows, columns int) mpc_core.RMat {
	out := make(mpc_core.RMat, rows)
	for i := range out {
		out[i] = values[i*columns : (i+1)*columns]
	}
	return out
}

// batchMatMul evaluates independent matrix products with one Beaver batch.
func (ast *AssocTest) batchMatMul(left, right []mpc_core.RMat) []mpc_core.RMat {
	if len(left) != len(right) {
		panic("batched matrix count mismatch")
	}
	rtype := ast.general.mpcObj[0].GetRType()
	out := make([]mpc_core.RMat, len(left))
	rows := make([]int, len(left))
	inner := make([]int, len(left))
	columns := make([]int, len(left))
	var leftFlat, rightFlat mpc_core.RVec
	for gene := range left {
		rightRows := 0
		rows[gene], inner[gene] = matrixDims(left[gene])
		rightRows, columns[gene] = matrixDims(right[gene])
		if rows[gene] == 0 || rightRows == 0 || inner[gene] == 0 || columns[gene] == 0 {
			continue
		}
		if inner[gene] != rightRows {
			panic("batched matrix dimensions mismatch")
		}
		out[gene] = mpc_core.InitRMat(rtype.Zero(), rows[gene], columns[gene])
		for i := range left[gene] {
			leftFlat = append(leftFlat, left[gene][i]...)
		}
		for i := range right[gene] {
			rightFlat = append(rightFlat, right[gene][i]...)
		}
	}
	if len(leftFlat) == 0 {
		return out
	}
	mpcObj := ast.general.mpcObj[0]
	leftRevealed, leftMask := mpcObj.BeaverPartitionVec(leftFlat)
	rightRevealed, rightMask := mpcObj.BeaverPartitionVec(rightFlat)
	leftOffset, rightOffset := 0, 0
	var products mpc_core.RVec
	for gene := range out {
		if out[gene] == nil {
			continue
		}
		leftSize := rows[gene] * inner[gene]
		rightSize := inner[gene] * columns[gene]
		product := mpcObj.BeaverMultMat(
			matrixFromFlat(leftRevealed[leftOffset:leftOffset+leftSize], rows[gene], inner[gene]),
			matrixFromFlat(leftMask[leftOffset:leftOffset+leftSize], rows[gene], inner[gene]),
			matrixFromFlat(rightRevealed[rightOffset:rightOffset+rightSize], inner[gene], columns[gene]),
			matrixFromFlat(rightMask[rightOffset:rightOffset+rightSize], inner[gene], columns[gene]),
		)
		for i := range product {
			products = append(products, product[i]...)
		}
		leftOffset += leftSize
		rightOffset += rightSize
	}
	products = mpcObj.BeaverReconstructVec(products)
	products = mpcObj.TruncVec(products, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	index := 0
	for gene := range out {
		for i := range out[gene] {
			copy(out[gene][i], products[index:index+len(out[gene][i])])
			index += len(out[gene][i])
		}
	}
	return out
}

func (ast *AssocTest) batchElemMul(left, right []mpc_core.RMat) []mpc_core.RMat {
	if len(left) != len(right) {
		panic("batched matrix count mismatch")
	}
	rtype := ast.general.mpcObj[0].GetRType()
	out := make([]mpc_core.RMat, len(left))
	var a, b mpc_core.RVec
	for gene := range left {
		rows, columns := matrixDims(left[gene])
		rightRows, rightColumns := matrixDims(right[gene])
		if rows != rightRows || columns != rightColumns {
			panic("batched element dimensions mismatch")
		}
		out[gene] = mpc_core.InitRMat(rtype.Zero(), rows, columns)
		for i := 0; i < rows; i++ {
			a = append(a, left[gene][i]...)
			b = append(b, right[gene][i]...)
		}
	}
	if len(a) == 0 {
		return out
	}
	products := ast.ssMul(a, b)
	index := 0
	for gene := range out {
		for i := range out[gene] {
			copy(out[gene][i], products[index:index+len(out[gene][i])])
			index += len(out[gene][i])
		}
	}
	return out
}

func (ast *AssocTest) batchScaleRows(rows []mpc_core.RVec, matrices []mpc_core.RMat) []mpc_core.RMat {
	if len(rows) != len(matrices) {
		panic("batched row scale count mismatch")
	}
	rtype := ast.general.mpcObj[0].GetRType()
	scales := make([]mpc_core.RMat, len(rows))
	for gene := range rows {
		n, columns := matrixDims(matrices[gene])
		if n != len(rows[gene]) {
			panic("batched row scale dimensions mismatch")
		}
		scales[gene] = mpc_core.InitRMat(rtype.Zero(), n, columns)
		for i := 0; i < n; i++ {
			for j := 0; j < columns; j++ {
				scales[gene][i][j] = rows[gene][i]
			}
		}
	}
	return ast.batchElemMul(scales, matrices)
}

func (ast *AssocTest) batchScaleMatrices(matrices []mpc_core.RMat, factor float64) []mpc_core.RMat {
	mpcObj := ast.general.mpcObj[0]
	out := make([]mpc_core.RMat, len(matrices))
	var flat mpc_core.RVec
	for gene, matrix := range matrices {
		rows, columns := matrixDims(matrix)
		out[gene] = mpc_core.InitRMat(mpcObj.GetRType().Zero(), rows, columns)
		for i := range matrix {
			flat = append(flat, matrix[i]...)
		}
	}
	if len(flat) == 0 {
		return out
	}
	flat = ast.ssPMul(flat, factor)
	index := 0
	for gene := range out {
		for i := range out[gene] {
			copy(out[gene][i], flat[index:index+len(out[gene][i])])
			index += len(out[gene][i])
		}
	}
	return out
}

func (ast *AssocTest) batchVectorDots(left, right []mpc_core.RVec) mpc_core.RVec {
	if len(left) != len(right) {
		panic("batched vector count mismatch")
	}
	rtype := ast.general.mpcObj[0].GetRType()
	out := mpc_core.InitRVec(rtype.Zero(), len(left))
	var a, b mpc_core.RVec
	for gene := range left {
		if len(left[gene]) != len(right[gene]) {
			panic("batched vector dimensions mismatch")
		}
		a = append(a, left[gene]...)
		b = append(b, right[gene]...)
	}
	if len(a) == 0 {
		return out
	}
	products := ast.ssMul(a, b)
	index := 0
	for gene := range left {
		for range left[gene] {
			out[gene] = out[gene].Add(products[index])
			index++
		}
	}
	return out
}

func (ast *AssocTest) batchMatrixDots(left, right []mpc_core.RMat) mpc_core.RVec {
	a := make([]mpc_core.RVec, len(left))
	b := make([]mpc_core.RVec, len(right))
	for gene := range left {
		rows, columns := matrixDims(left[gene])
		rightRows, rightColumns := matrixDims(right[gene])
		if rows != rightRows || columns != rightColumns {
			panic("batched matrix dot dimensions mismatch")
		}
		for i := 0; i < rows; i++ {
			a[gene] = append(a[gene], left[gene][i]...)
			b[gene] = append(b[gene], right[gene][i]...)
		}
	}
	return ast.batchVectorDots(a, b)
}

func (ast *AssocTest) batchRowDots(left, right []mpc_core.RMat) []mpc_core.RVec {
	products := ast.batchElemMul(left, right)
	rtype := ast.general.mpcObj[0].GetRType()
	out := make([]mpc_core.RVec, len(products))
	for gene, matrix := range products {
		out[gene] = mpc_core.InitRVec(rtype.Zero(), len(matrix))
		for i := range matrix {
			for _, value := range matrix[i] {
				out[gene][i] = out[gene][i].Add(value)
			}
		}
	}
	return out
}

func (ast *AssocTest) batchTraceProducts(left, right []mpc_core.RMat) mpc_core.RVec {
	a := make([]mpc_core.RVec, len(left))
	b := make([]mpc_core.RVec, len(right))
	for gene := range left {
		rows, columns := matrixDims(left[gene])
		rightRows, rightColumns := matrixDims(right[gene])
		if rows != rightColumns || columns != rightRows {
			panic("batched trace dimensions mismatch")
		}
		for i := 0; i < rows; i++ {
			for j := 0; j < columns; j++ {
				a[gene] = append(a[gene], left[gene][i][j])
				b[gene] = append(b[gene], right[gene][j][i])
			}
		}
	}
	return ast.batchVectorDots(a, b)
}

func traceShares(rtype mpc_core.RElem, matrices []mpc_core.RMat) mpc_core.RVec {
	out := mpc_core.InitRVec(rtype.Zero(), len(matrices))
	for gene, matrix := range matrices {
		for i := range matrix {
			out[gene] = out[gene].Add(matrix[i][i])
		}
	}
	return out
}
