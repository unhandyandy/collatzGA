
using StatsBase,Keras

import Keras.Layers: Dense, Activation

model = Sequential()
add!(model, Dense(80, input_dim=735))
add!(model, Activation(:relu))
add!(model, Dense(10))
add!(model, Activation(:softmax))

compile!(model; loss=:categorical_crossentropy, optimizer=:sgd, metrics=[:accuracy])

h = fit!(model, rand(1000, 735), rand(1000, 10); nb_epoch=100, batch_size=32, verbose=1)

evaluate(model, rand(10, 735), rand(10, 10); batch_size=5, verbose=1)

predict(model, rand(10, 735); batch_size=5, verbose=1)

