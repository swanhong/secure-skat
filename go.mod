module github.com/hhcho/sfgwas

go 1.23.1

require (
	github.com/BurntSushi/toml v1.6.0
	github.com/aead/chacha20 v0.0.0-20180709150244-8b13a72661da
	github.com/hhcho/frand v0.0.0-00010101000000-000000000000
	github.com/hhcho/mpc-core v0.0.0-20240903135357-56a83f968a6a
	github.com/raulk/go-watchdog v1.3.0
	github.com/tuneinsight/lattigo/v6 v6.1.1
	go.dedis.ch/onet/v3 v3.2.10
	gonum.org/v1/gonum v0.9.3
)

require (
	github.com/ALTree/bigfloat v0.0.0-20220102081255-38c8b72a9924 // indirect
	github.com/benbjohnson/clock v1.3.0 // indirect
	github.com/containerd/cgroups v0.0.0-20201119153540-4cbc285b3327 // indirect
	github.com/coreos/go-systemd/v22 v22.1.0 // indirect
	github.com/davecgh/go-spew v1.1.1 // indirect
	github.com/daviddengcn/go-colortext v0.0.0-20180409174941-186a3d44e920 // indirect
	github.com/docker/go-units v0.4.0 // indirect
	github.com/elastic/gosigar v0.12.0 // indirect
	github.com/godbus/dbus/v5 v5.0.3 // indirect
	github.com/gogo/protobuf v1.3.1 // indirect
	github.com/google/go-cmp v0.5.8 // indirect
	github.com/kr/text v0.2.0 // indirect
	github.com/opencontainers/runtime-spec v1.0.2 // indirect
	github.com/pkg/errors v0.9.1 // indirect
	github.com/pmezard/go-difflib v1.0.0 // indirect
	github.com/stretchr/testify v1.8.0 // indirect
	golang.org/x/crypto v0.31.0 // indirect
	golang.org/x/exp v0.0.0-20230321023759-10a507213a29 // indirect
	golang.org/x/sys v0.28.0 // indirect
	golang.org/x/xerrors v0.0.0-20191204190536-9bdfabe68543 // indirect
	gopkg.in/yaml.v3 v3.0.1 // indirect
)

replace github.com/tuneinsight/lattigo/v6 => ../matrix_ckks

replace github.com/hhcho/frand => lukechampine.com/frand v1.5.1
