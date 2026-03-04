#!/bin/bash

TEST_NAME=${1:-TestSecureSKATEndToEnd}
NUM_MAIN_PARTY=2
LOG_PREFIX=test_stdout_party

pkill -f "skat_test.test" 2>/dev/null

echo "Compiling test binary..."
go test -c ./gwas/ -o skat_test.test || { echo "Compilation failed"; exit 1; }

rm -f ${LOG_PREFIX}*.txt

for (( i = 0; i <= $NUM_MAIN_PARTY; i++ ))
do
  echo "Running PID=$i"
  CMD="PID=$i ./skat_test.test -test.v -test.run $TEST_NAME > ${LOG_PREFIX}${i}.txt 2>&1"
  if [ $i = $NUM_MAIN_PARTY ]; then
    eval $CMD
  else
    eval $CMD &
    sleep 2
  fi
done

wait
echo "Tests completed. Check ${LOG_PREFIX}*.txt for outputs."
