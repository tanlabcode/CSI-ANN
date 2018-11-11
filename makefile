s2: CSIANNTrain.o OneLayersTimeDelayNeuralNet.o NeuralNetLib.o mem_loc.o
	gcc -o s2 CSIANNTrain.o OneLayersTimeDelayNeuralNet.o NeuralNetLib.o mem_loc.o -lm
CSIANNTrain.o: CSIANNTrain.c OneLayersTimeDelayNeuralNet.h
	gcc -c -g CSIANNTrain.c
OneLayersTimeDelayNeuralNet.o: OneLayersTimeDelayNeuralNet.c OneLayersTimeDelayNeuralNet.h
	gcc -c -g OneLayersTimeDelayNeuralNet.c
NeuralNetLib.o: NeuralNetLib.c NeuralNetLib.h
	gcc -c -g NeuralNetLib.c
mem_loc.o: mem_loc.c mem_loc.h
	gcc -c -g mem_loc.c
