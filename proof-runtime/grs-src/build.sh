F=/workspace/rh-research/tools/flintc/root/usr; L=/workspace/rh-research/tools/flintc/root/usr/lib/x86_64-linux-gnu
gcc $2 -fopenmp -I$F/include -I$F/include/x86_64-linux-gnu -o $1 $3 -L$L -lflint -lm -Wl,-rpath,$L
