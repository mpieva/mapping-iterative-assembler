set -e
./bootstrap.sh
./configure --prefix=/home/public/usr64
make
VERSION=`svnversion`
make install prefix=/home/public/usr64/stow/mia-$VERSION
cd /home/public/usr64/stow
stow -D mia-*
stow -v mia-$VERSION
cd -
