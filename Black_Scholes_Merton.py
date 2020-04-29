#Black-Scholes Model for Options Pricing
import io
import sys
import os
from scipy.stats import norm
import math
import numpy
import matplotlib.pyplot as plt

def black_scholes(price, strike_price, volatility,Time_to_Expiration,interest_rate,dividen_yield):

	#constants. Straight_forward
	K = float(strike_price)
	x = float(price)
	sig = float(volatility)
	t = float(Time_to_Expiration)
	r = float(interest_rate)
	d = float(dividen_yield)
	
	
	d1 = (1/(sig*math.sqrt(t)))*(numpy.log(x/K) + (r + (sig*sig)/2)*t)
	d2 = (1/(sig*math.sqrt(t)))*(numpy.log(x/K) + (r - (sig*sig)/2)*t)
	
	BSMC = x*math.exp(-d*t)*norm.cdf(d1) - K*math.exp(-r*t)*norm.cdf(d2)
	BSMP = K*math.exp(-r*t)*norm.cdf(-d2) - x*math.exp(-d*t)*norm.cdf(-d1)
	print("CALL: " ,BSMC)
	print("PUT: " ,BSMP)
	deltacall = norm.cdf(d1)
	
	#print("Delta Call: " , deltacall)
	
	Intrinsic = max(x -  K,0)
	Extrinsic_CALL=  BSMC - Intrinsic
	print("Extrinsic Value CALL: ",  Extrinsic_CALL)
	
	#With the value of the Extrinsic Value
	#we can calculate the effective Interest rate you pay for the remaining leverage
	#an effective calculation should take into account the dividends you miss my being in a call
	#rather than in shares
	
	Effective_Interest_Rate = Extrinsic_CALL/(K)
	print("Effective_Interest_Rate: ",Effective_Interest_Rate)
	
	#print("Delta Put: ", deltaput)
	return Effective_Interest_Rate;

def opLev(mu,sig,price,rate,div,time):
        m = float(mu)
        s = float(sig)
        p = float(price)
        r = float(rate)
        t = float(time)
        d = float(div)
        optimal_k = 0
        for i in range(1,int(p)):
                Leverage = p/(p-i)
                
                #text_trap = io.StringIO()
                #sys.stdout = text_trap
                x = black_scholes(p,i,s,t,r,d)
                #sys.stdout =  sys.__stdout__
                
                f = (m-x)/(s*s)
                if (f - Leverage <= 0.02):
                    optimal_k = i
                    break
        array = [optimal_k,black_scholes(p,optimal_k,s,t,r,d)]
        return array

def loop(price,vol,rate,div):
        #lo_p = max(int(price - price*vol*3),1)
        lo_p = 1
        hi_p = int(price + price*vol*3)
        p = float(price)
        v = float(vol)
        r_f = float(rate)
        r_div = float(div)
        x = []
        k = []
        r = []
        EL = []
        kelly = []
        
        for i in range(lo_p,hi_p,1):
                
                #print(i)
                text_trap = io.StringIO()
                sys.stdout = text_trap
                Interest = black_scholes(p,i,v,1,r_f,r_div) #type error: can only cocatenate list(not "float:) to list
                sys.stdout =  sys.__stdout__
                #print(Interest)
                x.append(Interest)
                k.append(i)
                r.append(0.017)
                kelly.append((0.07 - (Interest))/(0.18*0.18))
                if i < price:
                        EL.append(price/(price - i))
                elif i >= price:
                        EL.append(1) 
        plt.figure(1)
        
        #plt.subplot(211)
        #plt.scatter(k,x)
        #plt.scatter(k,r)
        ax = plt.subplot(211)
        ax.plot(k,x, label = "Effective Option Rate")
        ax.plot(k,r,label = "Risk-Free Rate")
        plt.axvline(x = price + v*price,color = "red",dashes = (1,2) )
        plt.axvline(x = price - v*price,color = "red",dashes = (1,2))
        plt.axvline(x = (price - v*2*price), color = "pink",dashes = (1,2))
        plt.axvline(x = (price + v*2*price),color = "pink",dashes = (1,2) )
        ax.legend()           
        ax2 = plt.subplot(212)
        ax2.plot(k,EL,label = "Option Leverage")
        ax2.plot(k,kelly,label = "Kelly Leverage")
        ax2.legend()
        plt.show()

        return ;
def monte_carlo(mean,vol,years):
        #inputs in yearly numbers
        mu = mean
        sig = vol
        #kelly = mu/(sig*sig) # second order approximation
        #print(kelly)
        kp = []
        x0 = 1
        x1 = 1
        price = []
        price.append(x0)
        p2 = []
        p2.append(x1)
        kp.append(x1)
        dt = 1/252
        total_steps = years*252
        y = []
        y.append(0.8)
        ks = 0.8
        #Geometric Brownian Motion
                
        for i in range(total_steps):
                y.append(0.8)
                if i%(252*2) == 0 :
                        x0 = 1
                        ks = 0.8*kp[i]
                
                jump = numpy.exp(numpy.random.poisson(0.002)*-0.05)
                pc = numpy.exp(dt*(mu - sig*sig/2))*numpy.exp(numpy.random.normal(0,dt**(1/2))*sig)*jump
                
                x1 = x1*pc
                if x0 > 0.8 :
                        
                        ky = kp[i]  +  (x1 - p2[i])
                else:
                         ky = ks
                x0 = x0*pc
                kp.append(ky)
                price.append(x0)
                p2.append(x1)

        t = numpy.linspace(0,total_steps, total_steps + 1)
        plt.figure(1)
        plt.subplot(211)
        plt.plot(t,price, color = "blue")
        plt.plot(y, color = "black")
        #plt.subplot(212)
        #plt.plot(t, kp, color = "black")
        for i in range(years + 1):
                #plt.subplot(211)
                plt.axvline(x = 252*2*i,color = "red",dashes = (1,2) )
        plt.subplot(212)
        plt.plot(t,p2,color = "blue")
        plt.plot(t,kp,color = "black")
        for i in range(years + 1):
                #plt.subplot(211)
                plt.axvline(x = 252*2*i,color = "red",dashes = (1,2) )
        plt.show()

        return
        
