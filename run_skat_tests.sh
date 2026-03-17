#!/bin/bash

TEST_NAME=${1:-TestSecureSKATEndToEnd}
RUN_SUFFIX=${2:-}
NUM_MAIN_PARTY=2
LOG_PREFIX=test_stdout

if [ -n "$RUN_SUFFIX" ]; then
  LOG_PREFIX="${LOG_PREFIX}_${RUN_SUFFIX}_party"
else
  LOG_PREFIX="${LOG_PREFIX}_party"
fi

pkill -f "skat_test.test" 2>/dev/null

echo "Compiling test binary..."
go test -c ./gwas/ -o skat_test.test || { echo "Compilation failed"; exit 1; }

rm -f ${LOG_PREFIX}*.txt

for (( i = 0; i <= $NUM_MAIN_PARTY; i++ ))
do
  echo "Running PID=$i"
  CMD="PID=$i TEST_RUN_SUFFIX=$RUN_SUFFIX ./skat_test.test -test.v -test.run $TEST_NAME > ${LOG_PREFIX}${i}.txt 2>&1"
  if [ $i = $NUM_MAIN_PARTY ]; then
    eval $CMD
  else
    eval $CMD &
    sleep 2
  fi
done

wait
if [ -n "$RUN_SUFFIX" ]; then
  echo "Tests completed. Logs: ${LOG_PREFIX}*.txt"
  echo "Outputs: out/party*_${RUN_SUFFIX}, cache/party*_${RUN_SUFFIX}"
else
  echo "Tests completed. Check ${LOG_PREFIX}*.txt for outputs."
fi
