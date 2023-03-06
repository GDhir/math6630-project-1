import numpy as np
import matplotlib.pyplot as plt
import math

def computeuexact( x ):

    return math.sin( math.pi*x )

def computek( x ):

    kval = 2

    for l in range( 1, 6):

        kval += ( math.sin( l*math.pi*x ) )/( l + 1 )

    return kval

def computef( x ):

    fval = -2*( (math.pi)**2 )*math.sin( math.pi*x )

    for l in range(1, 6):

        fval = fval - ( (math.pi)**2 )*math.sin( ( math.pi ) * x )*( ( math.sin( l*math.pi*x )/( l + 1 ) ) ) \
        + l*( math.pi**2 )*math.cos( math.pi*x )*( math.cos( l*math.pi*x )/( l + 1 ) )

    return fval


def findsol( N ):

    h = 1/( N - 1 )

    error = 1e5

    tol = 1e-6

    u = np.zeros( N )
    uprev = np.zeros( N )

    while error > tol:

        error = 0

        for j in range( 1, N - 1 ):

            xj = j*h

            kj1 = computek( xj + h/2 )
            kj2 = computek( xj - h/2 )

            A = kj1 + kj2
            fj = computef( xj )

            u[ j ] = ( ( kj1*uprev[ j + 1 ] + kj2*u[ j - 1 ] )/( h**2 ) - fj )*(h**2)/A

            # error += ( uprev[j] - u[j] )**2

            uprev[j] = u[j]

        for j in range( 1, N - 1 ):

            xj = j*h
            fj = computef( xj )
            kj1 = computek( xj + h/2 )
            kj2 = computek( xj - h/2 )
            error += ( ( kj1*u[ j + 1 ] + kj2*u[ j - 1 ] - ( kj1 + kj2 )*u[j] )/( h**2 ) - fj )**2

        error = math.sqrt( error )/N
        # print(error)

    return u

def runscript():

    Nvals = [ 17, 25, 33, 49, 65 ]
    hvals = []

    errorvals = np.zeros( len( Nvals ) )

    for idx, N in enumerate(Nvals):

        h = 1/( N - 1 )
        u = findsol(N)
        x = [0]

        for j in range( 1, N - 1 ):

            xj = j*h
            x.append( xj )
            uexact = computeuexact( xj )
            errorvals[idx] += ( uexact - u[j] )**2

        errorvals[idx] = np.sqrt( errorvals[idx]/N )

        x.append(1)

        hvals.append( h**2 )

        if N == 65:

            plt.figure()
            plt.plot( x, u, "-o" )
            plt.xlabel( "X" )
            plt.ylabel( "Solution (u)" )
            plt.savefig( "Solution_Plot_65.png" )

    slope_intercept = np.polyfit( np.log2(Nvals), np.log2(errorvals), 1 )

    print(slope_intercept)




    plt.figure()
    plt.loglog( Nvals, errorvals, "-o", label="Solution Error" )
    plt.loglog( Nvals, hvals, "-o", label = "$h^2$ Plot" )
    plt.xlabel( "Number of Points in the Grid" )
    plt.ylabel( "Error or $h^2$" )
    plt.legend()
    plt.savefig( "Error_Plot.png" )
    plt.show()




if __name__ == "__main__":
    
    # u = findsol( 65 )

    # plt.plot( u )

    # plt.show()

    runscript()