git clone https://github.com/cran/rda.git
R CMD INSTALL --no-multiarch --with-keep.source rda
rm -rd rda
git clone https://github.com/longhaiSK/HTLR.git
cd HTLR
git checkout legacy
cd ..
R CMD INSTALL --no-multiarch --with-keep.source HTLR
rm -rf HTLR