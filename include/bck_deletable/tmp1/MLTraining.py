import numpy as np
import sys
import time as ct
import os.path
import math
import torch
import torch.nn as nn
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
from sklearn.preprocessing import MinMaxScaler

def DataNormlize(model, inp, grad, tangent):
    inp_scaler = MinMaxScaler()
    inp_scaled = inp_scaler.fit_transform(inp)
    scales_inp = inp_scaler.scale_
    limits_inp = inp_scaler.min_
    X_train = inp_scaled[0:2200,:].astype(np.float32)
    X_test = inp_scaled[2200:,:].astype(np.float32)

    grad_scaler = MinMaxScaler()
    grad_scaled = grad_scaler.fit_transform(grad) 
    scales_grad = grad_scaler.scale_
    limits_grad =  grad_scaler.min_
    y_train = grad_scaled[0:2200,:].astype(np.float32)
    y_test = grad_scaled[2200:,:].astype(np.float32)

    tangent_scaler = MinMaxScaler()
    tangent_scaled = tangent_scaler.fit_transform(tangent) 
    scales_tangent = tangent_scaler.scale_
    limits_tangent =  tangent_scaler.min_
    tangent_train = tangent_scaled[0:2200,:].astype(np.float32)
    tangent_test = tangent_scaled[2200:,:].astype(np.float32)
    
    print(inp_scaled)
    #print((inp-limits_inp)/scales_inp)
    #print()
    #grad_scaled = grad*scales_grad + limits_grad
    print("limits_inp\n",limits_inp)
    print("scales_inp\n",scales_inp)
    print("Normalize Example")
    print("input (ev, es): 0.2, 0.3")
    print("scaled_input = inp*scales_inp+limits_inp")
    print("inp1 = 0.2*",scales_inp[0],"+",limits_inp[0])
    print("inp2 = 0.3*",scales_inp[1],"+",limits_inp[1])
    model.inp_scaler = inp_scaler
    model.limits_inp = limits_inp
    model.scales_inp = scales_inp
    model.grad_scaler = grad_scaler
    model.scales_grad = scales_grad
    model.limits_grad = limits_grad

    X_train = torch.from_numpy(X_train)
    X_train.requires_grad= True
    X_test = torch.from_numpy(X_test)
    X_test.requires_grad = True
    y_train = torch.from_numpy(y_train)
    y_test = torch.from_numpy(y_test)    
    return X_train, X_test, y_train, y_test

def Training(model, E, nu):
    DataGen(E, nu)
    inp, grad, tangent = ReadData()  
    print(inp)
    X_train, X_test, y_train, y_test = DataNormlize(model, inp, grad, tangent)
    ModelTraining(model, 0.0005, 1000, X_train, y_train, X_test, y_test)


class NeuralNetwork(nn.Module):
    
    ##inp_scaler
    limits_inp: np.array
    scales_inp: np.array
    #grad_scaler
    scales_grad: np.array
    limits_grad: np.array
        
    def __init__(self, NetSize=[2,64,128,2]):
        super(NeuralNetwork,self).__init__()
        self.d1 = nn.Linear(NetSize[0],NetSize[1])
        self.act1 = nn.ReLU()
        self.d2 = nn.Linear(NetSize[1],NetSize[2])
        self.act2 = nn.ReLU()
        self.d3 = nn.Linear(NetSize[2],NetSize[3])
        
        
    def forward(self,x):
        x = self.d1(x)
        x = self.act1(x)
        x = torch.mul(x,x)
        x = torch.mul(x,x)
        x = self.d2(x)
        x = self.act2(x)
        x = torch.mul(x,x)
        out = self.d3(x)
        return out
    
    
    
def DataGen(E, nu):
    # data range of ev, es
    min_ev = -2.00
    max_ev = 2.00
    min_es = 0. # es is always positive
    max_es = 3.00

    num_points = 50

    # define material parameters
    #E = E
    #nu = nu
    K_bulk = E/3. / (1. - 2. * nu)
    G_shear = E/ 2./ (1. + nu)

    Xgrid, Ygrid = np.meshgrid(np.linspace(min_ev, max_ev, num_points), np.linspace(min_es, max_es, num_points))
    XX = Xgrid.reshape(-1,1) # X coordinate point set for volumetric strain
    YY = Ygrid.reshape(-1,1) # Y coordinate point set for dev strain
    np.random.shuffle(XX)
    np.random.shuffle(YY)

    def getf(X,Y):
        """
        compute the elastic storage energy in invariant forms, inputs:
            X is the volumetric strain
            Y is the deviatoric strain
        returns:
        the elastic function value psi
        """
        return 0.5 * K_bulk * X * X   +  G_shear * Y * Y * 3./2.

    def getfp(X,Y):
        p = K_bulk * X
        q = 3*G_shear * Y
        return p, q

    psi = getf(XX,YY)
    p,q = getfp(XX, YY)
    C11 = np.ones([p[:,0].size,1])*K_bulk
    C22 = np.ones([p[:,0].size,1])*3*G_shear
    # save all data
    dataset = pd.DataFrame({'eps11': XX[:,0], 'eps22': YY[:,0], 'psi_e': psi[:,0], 'W_e11': p[:,0], 'W_e22': q[:,0], 'C11': C11[:,0], 'C22': C22[:,0]})
    dataset.to_csv("hyperelastic_data.csv",index = False)
    return


def ReadData():
    DataFrame = pd.read_csv('hyperelastic_data.csv')
    # ------------------------------
    # Scale the data to range(0,1)
    # ------------------------------
    # scale inputs
    inp = DataFrame[["eps11", "eps22"]]


    # scale output
    out = DataFrame[["psi_e"]]
    grad = DataFrame[["W_e11","W_e22"]].values.reshape(-1,2)
    tangent = DataFrame[["C11","C22"]].values.reshape(-1,2)
    
    return inp, grad, tangent



def ModelTraining(model, lr, num_epochs, X_train, y_train, X_test, y_test):
    # train the model with the data
    loss_history = []
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    for epoch in range(num_epochs):

        # initialize the loss
        train_running_loss = 0.0
        test_loss = 0.0
        
        # compute y_pred and dX_pred = (d y_pred)/(d X_train)
        y_pred = model(X_train)
        external_grad = torch.ones_like(y_pred)
        y_pred.backward(gradient=external_grad, retain_graph=True, create_graph=True)
        dX_pred = X_train.grad 
        
        # compute the loss from the predictions and the ground truth
        loss = criterion(y_pred, y_train)
        #loss = criterion(dX_pred, dX_train) + criterion(y_pred, y_train)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        train_running_loss = loss.detach().item()

        # compute the test loss in the same way
        y_pred = model(X_test)
        external_grad = torch.ones_like(y_pred)
        #print(external_grad )
        #exit(1)
        y_pred.backward(gradient=external_grad, retain_graph=False)
        dX_pred = X_test.grad
        test_loss = criterion(y_pred, y_test)
        X_test.grad.data.zero_()
        
        #loss_history.append([loss, test_loss])
        loss_history.append([loss.detach().numpy(),test_loss.detach().numpy()])
        # print the loss every 100 epochs
        if epoch%100==0:    
            print('Epoch: %d | Loss = %.8f | Test Loss = %.8f'%(epoch, train_running_loss, test_loss))
    #figure        
    plt.plot(loss_history)

def PlotTrainingResult(model, E, nu):
    K_bulk = E/3. / (1. - 2. * nu)
    G_shear = E/ 2./ (1. + nu)

    min_ev = -0.20
    max_ev = 0.20
    min_es = 0.
    max_es = 0.30
    num_test_points = 60

    X_plot, Y_plot = np.meshgrid(np.linspace(min_ev, max_ev, num_test_points), np.linspace(min_es, max_es, num_test_points))
    test_inp = np.concatenate((X_plot.reshape(-1,1), Y_plot.reshape(-1,1)), axis=1)
    p_true = K_bulk * test_inp[:,0]
    q_true = 3*G_shear * test_inp[:,1]


    fig, ax1 = plt.subplots(1, 2, figsize=(8, 3))
    ax = ax1[0]
    c = ax.pcolor(X_plot, Y_plot, p_true.reshape(num_test_points,num_test_points), cmap='jet')
    ax.set_title('p true', fontsize=16)
    fig.colorbar(c, ax=ax)
    #plt.show()


    grad_pred = model.forward(torch.from_numpy(model.inp_scaler.transform(test_inp.astype(np.float32))))

    #print(p_pred.detach().numpy().size)
    grad_pred = grad_pred.detach().numpy()
    grad_pred = (grad_pred - model.limits_grad)/model.scales_grad
    #print(size(p_pred.detach().numpy()))

    #print(model.forward(torch.from_numpy(test_inp[0,:].astype(np.float32))))
    ax = ax1[1]
    c = ax.pcolor(X_plot, Y_plot, grad_pred[:,0].reshape(num_test_points,num_test_points), cmap='jet')
    ax.set_title('p pred', fontsize=16)
    fig.colorbar(c, ax=ax)
    plt.show()
    
    fig, ax1 = plt.subplots(1, 2, figsize=(8, 3))
    ax = ax1[0]
    c = ax.pcolor(X_plot, Y_plot, q_true.reshape(num_test_points,num_test_points), cmap='jet')
    ax.set_title('q true', fontsize=16)
    fig.colorbar(c, ax=ax)
    #plt.show()


    grad_pred = model.forward(torch.from_numpy(model.inp_scaler.transform(test_inp.astype(np.float32))))

    #print(p_pred.detach().numpy().size)
    grad_pred = grad_pred.detach().numpy()
    grad_pred = (grad_pred - model.limits_grad)/model.scales_grad
    #print(size(p_pred.detach().numpy()))

    #print(model.forward(torch.from_numpy(test_inp[0,:].astype(np.float32))))
    ax = ax1[1]
    c = ax.pcolor(X_plot, Y_plot, grad_pred[:,1].reshape(num_test_points,num_test_points), cmap='jet')
    ax.set_title('q pred', fontsize=16)
    fig.colorbar(c, ax=ax)
    plt.show()
