package crypto

import (
	"bufio"
	"encoding/binary"
	"io"
	"os"
)

func MarshalCM(cm CipherMatrix) ([]byte, []byte) {
	cmBytes, ctSizes, err := cm.MarshalBinary()
	if err != nil {
		panic(err)
	}

	r, c := len(cm), len(cm[0])
	intsize := uint64(8)
	sizesbuf := make([]byte, intsize*uint64(r)*uint64(c))

	offset := uint64(0)
	for i := range ctSizes {
		for j := range ctSizes[i] {
			binary.LittleEndian.PutUint64(sizesbuf[offset:offset+intsize], uint64(ctSizes[i][j]))
			offset += intsize
		}
	}

	return sizesbuf, cmBytes
}

func UnmarshalCM(cryptoParams *CryptoParams, r, c int, sbytes, ctbytes []byte) CipherMatrix {
	intsize := uint64(8)
	offset := uint64(0)
	sizes := make([][]int, r)
	for i := range sizes {
		sizes[i] = make([]int, c)
		for j := range sizes[i] {
			sizes[i][j] = int(binary.LittleEndian.Uint64(sbytes[offset:]))
			offset += intsize
		}
	}

	var cm CipherMatrix
	if err := (&cm).UnmarshalBinary(cryptoParams, ctbytes, sizes); err != nil {
		panic(err)
	}

	return cm
}

func SaveCipherMatrixToFile(_ *CryptoParams, cm CipherMatrix, filename string) {
	file, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	sbytes, cmbytes := MarshalCM(cm)

	nrbuf := make([]byte, 4)
	ncbuf := make([]byte, 4)
	binary.LittleEndian.PutUint32(nrbuf, uint32(len(cm)))
	binary.LittleEndian.PutUint32(ncbuf, uint32(len(cm[0])))

	sbuf := make([]byte, 8)
	cmbuf := make([]byte, 8)
	binary.LittleEndian.PutUint64(sbuf, uint64(len(sbytes)))
	binary.LittleEndian.PutUint64(cmbuf, uint64(len(cmbytes)))

	if _, err := writer.Write(nrbuf); err != nil {
		panic(err)
	}
	if _, err := writer.Write(ncbuf); err != nil {
		panic(err)
	}
	if _, err := writer.Write(sbuf); err != nil {
		panic(err)
	}
	if _, err := writer.Write(sbytes); err != nil {
		panic(err)
	}
	if _, err := writer.Write(cmbuf); err != nil {
		panic(err)
	}
	if _, err := writer.Write(cmbytes); err != nil {
		panic(err)
	}

	if err := writer.Flush(); err != nil {
		panic(err)
	}
}

func LoadCipherMatrixFromFile(cps *CryptoParams, filename string) CipherMatrix {
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	reader := bufio.NewReader(file)
	ibuf := make([]byte, 4)
	if _, err := io.ReadFull(reader, ibuf); err != nil {
		panic(err)
	}
	nrows := int(binary.LittleEndian.Uint32(ibuf))
	if _, err := io.ReadFull(reader, ibuf); err != nil {
		panic(err)
	}
	numCtxPerRow := int(binary.LittleEndian.Uint32(ibuf))

	sbuf := make([]byte, 8)
	if _, err := io.ReadFull(reader, sbuf); err != nil {
		panic(err)
	}
	sbyteSize := binary.LittleEndian.Uint64(sbuf)
	sdata := make([]byte, sbyteSize)
	if _, err := io.ReadFull(reader, sdata); err != nil {
		panic(err)
	}

	cmbuf := make([]byte, 8)
	if _, err := io.ReadFull(reader, cmbuf); err != nil {
		panic(err)
	}
	cbyteSize := binary.LittleEndian.Uint64(cmbuf)
	cdata := make([]byte, cbyteSize)
	if _, err := io.ReadFull(reader, cdata); err != nil {
		panic(err)
	}

	return UnmarshalCM(cps, nrows, numCtxPerRow, sdata, cdata)
}
