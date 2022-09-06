# Generate the figures needed for modecoupling.

# Make Obsolete the old ones.
for file in `fgrep includegraphics modecoupling.tex | sed -e 's/^[^{]*{//' -e 's/\}.*$//'`; do mv $file* obsolete/ ; done

# Delete old ones (now irrelevant)
for file in `fgrep includegraphics modecoupling.tex | sed -e 's/^[^{]*{//' -e 's/\}.*$//'`; do rm -f $file.pdf; done

make continmodez
make osolvomk
make ocontgen

./continmodez
mv plot0001.ps cmodez.ps
epstopdf cmodez.ps

#./osolvomk -oc1 -k.20 -no20 -nk4 -p.25 -M -w-3
#mv plot0001.ps F2004B1p025.ps
#epstopdf F2004B1p025.ps
#mv plot0002.ps F2004B1p025d.ps
#epstopdf F2004B1p025d.ps

# Higher resolution
./osolvomk -oc1 -k.20 -no40 -nk4 -p.25 -M -w-3
mv plot0001.ps F4004B1p025.ps
epstopdf F4004B1p025.ps
mv plot0002.ps F4004B1p025d.ps
epstopdf F4004B1p025d.ps

# Alternative to the above, not used:
#./osolvomk -no20 -nk6 -p.09 -k.3 -M -w-3
#mv plot0001.ps F2006B1p009.ps
#epstopdf F2006B1p009.ps
#mv plot0002.ps F2006B1p009d.ps
#bepstopdf F2006B1p009d.ps

mpiexec -n 8 ./ocontgen -oc.2 -p.25 -k.1 -c
mv plot0004.ps B02k1p25.ps
epstopdf B02k1p25.ps

mpiexec -n 8 ./ocontgen -oc.5 -p.25 -k.4 -c
mv plot0002.ps B05k4p25.ps
epstopdf B05k4p25.ps
mv plot0003.ps B05k4p25modes.ps
epstopdf B05k4p25modes.ps

mpiexec -n 8 ./ocontgen -oc1. -p.25 -k.1 -c
mv plot0002.ps B10k10p25.ps
epstopdf  B10k10p25.ps

cp modeshape.plt modeshape1628.plt
plottraces modeshape1628.plt totalmode1628.plt -pf-3
mv plot0001.ps modeshape1628.ps
epstopdf modeshape1628.ps

#High B plot
mpiexec -n 8 ./ocontgen -p.49 -k.15 -oc9  -or.05 -oi.001 -c
mv plot0004.ps B99k15p49.ps
epstopdf B99k15p49.ps

# Extern figure got from auxmodeversion (not obsoleted)
cp extern.ps extern3.ps
epstopdf extern3.ps

