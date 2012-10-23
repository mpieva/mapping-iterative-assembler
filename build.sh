set -e
./bootstrap.sh
./configure --prefix=/home/public/usr
make
VERSION=`git describe --always`
make install prefix=/home/public/usr/stow/mia-$VERSION
cd /home/public/usr/stow
stow -D mia-*
stow -v mia-$VERSION
cd -
