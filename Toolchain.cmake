#set(Boost_USE_STATIC_LIBS ON) #use this option if have BOOST static libraries (unusual)
set(BOOST_ROOT "/home/petra1/work/installs/boost_1_57_0/build/")

#in this particular case, libscalapack.a  includes blacs. You may need to specify libblacs.a separetely 
set(SCALAPACK_LIBRARIES "-L/home/petra1/work/installs/scalapack_installer/install/lib -lscalapack")