java -Dpams.home=./conf  -cp lib/epis.jar:lib/BehaviorSpace.jar:lib/jogl.jar:lib/MRJAdapter.jar:lib/jonas-client.jar:lib/NetLogo.jar:lib/jonas-cluster-daemon.jar:lib/Pams_properties.jar:./lib/jonas-generators-clientstubs.jar:./lib/asm-3.1.jar:./lib/jonas-generators-raconfig.jar:./lib/asm-commons-3.1.jar:./lib/log4j-1.2.14.jar:./lib/asm-util-3.1.jar:./lib/mysql-connector-java-5.1.6-bin.jar:./lib/client.jar:./lib/netlogoDriver.jar:./lib/picocontainer.jar:./lib/gluegen-rt.jar:./lib/quaqua.jar:./lib/jaxb-xjc.jar:./lib/scala-library.jar:./lib/jhotdraw.jar:./lib/smartclient.jar:./lib/jmf.jar:./lib/swing-layout.jar ummisco.cluster.loader.Load $1 $2

./pr $1 $2
