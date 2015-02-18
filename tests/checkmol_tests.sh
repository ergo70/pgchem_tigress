#!/bin/sh

PASS="1:T"
FAIL="1:F"
FAILS=0
RESULT=`matchmol -x molfiles/test2jbug.mol molfiles/test2jbug.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 1 passed';
else
	echo 'Test 1 failed';
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol -xs molfiles/test2jbug.mol molfiles/test2jbug.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 2 passed'
else
	echo 'Test 2 failed'
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol molfiles/test2jbug.mol molfiles/test2jbug.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 3 passed';
else
	echo 'Test 3 failed';
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol -s molfiles/test2jbug.mol molfiles/test2jbug.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 4 passed'
else
	echo 'Test 4 failed'
	FAILS=`expr ${FAILS} + 1`;
fi

if [ $FAILS -ne 0 ]; then
	echo 'Your Version of checkmol/matchmol contains a bug that was fixed in 0.2k! This affects equality detection, so an update is recommended.'
fi

FAILS=0
RESULT=`matchmol -x molfiles/compsan.mol molfiles/compsan.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 5 passed';
else
	echo 'Test 5 failed';
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol -xs molfiles/compsan.mol molfiles/compsan.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 6 passed'
else
	echo 'Test 6 failed'
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol molfiles/compsan.mol molfiles/compsan.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 7 passed';
else
	echo 'Test 7 failed';
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol -s molfiles/compsan.mol molfiles/compsan.mol`

if [ $RESULT = $PASS ]; then
	echo 'Test 8 passed'
else
	echo 'Test 8 failed'
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol molfiles/cyclooctane.mol molfiles/cis-pinonic-acid.mol`

if [ $RESULT = $FAIL ]; then
	echo 'Test 9 passed'
else
	echo 'Test 9 failed'
	FAILS=`expr ${FAILS} + 1`;
fi

RESULT=`matchmol molfiles/cycloheptane.mol molfiles/spirohexane-carboxylic-acid.mol`

if [ $RESULT = $FAIL ]; then
	echo 'Test 10 passed'
else
	echo 'Test 10 failed'
	FAILS=`expr ${FAILS} + 1`;
fi

if [ $FAILS -ne 0 ]; then
	echo 'Your Version of checkmol/matchmol does not compare correctly. Maybe the compiler produces wrong code. Try compiling without any optimization and try again'
fi

