from math import pi, e, log, sqrt, sin, cos

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

class UtilityCurve():
    """
    UtilityCurve is the Utility Curve representation that will be used.
    Although the underlying Utility Function can be called
    directly, it's preferred that value(deltaX) be called, and that when 
    arbitration is final that the progressByDeltaX(...) is used to advance
    the progression along the Utility Curve.
    """
    def __init__( self, x_0, x_25, x_50, x_75, x_100
                , calculate=True
                , weight=1.0
                , bounty=1.0
                ):
        """
        Initialize UtilityCurve.
        calculate may be set to false to bypass CDF generation,
        in which case x-points arguments are ignored and k-points are 
        not calculated.  They must be set before use, but this can
        be used to prevent recalculating UtilityCurveFunction at runtime.
        """
        self.true_min = x_0
        self.true_max = x_100
        x_0   = self.normalize(x_0)
        x_25  = self.normalize(x_25)
        x_50  = self.normalize(x_50)
        x_75  = self.normalize(x_75)
        x_100 = self.normalize(x_100)
        
        (   self.x_0 
          , self.x_50 
          , self.x_100 
          , self.k_1 
          , self.k_2 
          , self.k_3
          , self.k_4 ) = GenerateCDFConstructor(x_0, x_25, x_50, x_75, x_100)

        self.bounty = bounty
        self.weight = weight

        self.previousBaseMargin = False
        self.aip = False
        self.dmp = False
        self.inflectionP = False
        self.MSRcache = dict()


    def normalize(self, x):
        """
        Normalize value accoring to x_value init params
        """
        return float(x - self.true_min) / (self.true_max - self.true_min)

    def denormalize(self, x):
        """
        Denormalize value accoring to x_value init params
        """
        return float(x) * (self.true_max - self.true_min) + self.true_min
        
    def setBounty(self, bounty):
        self.bounty = bounty
        return self

    def setWeight(self, weight):
        self.weight = weight
        return self

    def __call__(self, x):
        """
        Returns value of UtilityCurve at normalized progression, or 
        affected by weight, and bounty, and denormalized.

        This is the easiest way to call the UtilityCurve.
        """
        return ( self.weight * self.bounty 
               * self.denormalize ( CDF( self.x_0, self.x_50, self.x_100
                                       , self.k_1, self.k_2, self.k_3, self.k_4
                                       , self.normalize(x)
               )                   )   )

                                      

    def raw_utility(self, x):
        """
        Returns Utility at x

        Not Normalized
        """
        return CDF(self.x_0, self.x_50, self.x_100, self.k_1, self.k_2, 
                self.k_3, self.k_4, x)

    def raw_margin(self, x):
        """
        Returns first derivative (Martin) at x

        Not Normalized
        """
        return CDF_prime(self.x_0, self.x_50, self.x_100, self.k_1, self.k_2, 
                self.k_3, self.k_4, x)

    def raw_growth(self, x):
        """
        Returns second derivative (Growth) at x

        Not Normalized
        """
        return CDF_prime_prime(self.x_0, self.x_50, self.x_100, self.k_1, 
                self.k_2, self.k_3, self.k_4, x)

    #----------------------------------------------------------------------
    #        For display now, efficient storage is not needed.
    #----------------------------------------------------------------------

    def getDict(self):
        return dict([ ('x_0'  , self.x_0),   ('x_50' , self.x_50)
                    , ('x_100', self.x_100), ('k_1'  , self.k_1 )
                    , ('k_2'  , self.k_2 ),  ('k_3'  , self.k_3 )
                    , ('k_4'  , self.k_4 )
                    ])

    def __repr__(self):
            return str(self.getDict())


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def CDF(x_0, x_50, x_100, k_1, k_2, k_3, k_4, x):
    """
    A Cumulative Distribution Function.
    Function is bound by [x_0, x_100]
    x_50 is the midpoint, k_1 through k_4 are transforming 
    variables.
    x will be a value on [0.0, 1.0]
    """
    if   x < x_0  : return 0.0
    elif x > x_100: return 1.0
    return k_3 / (1 + (k_1 * e ** (-k_2 * (x - x_50)))) + k_4

def CDF_prime(x_0, x_50, x_100, k_1, k_2, k_3, k_4, x):
    """
    the first derivative of a CDF function (margin)
    """
    if   x < x_0  : return 0.0
    elif x > x_100: return 0.0
    numerator =  k_1 * k_2 * k_3 * e ** (k_2 * (x - x_50))
    denominator = (k_1 + e ** (k_2 * (x - x_50))) ** 2
    return numerator / denominator

def CDF_prime_prime(x_0, x_50, x_100, k_1, k_2, k_3, k_4, x):
    """
    the second derivative of a CDF function (growth)
    """
    if   x < x_0  : return 0.0
    elif x > x_100: return 0.0
    num_0 = 2 * k_1 ** 2 * k_2 ** 2 * e ** (-2 * k_2 * (x - x_50))
    den_0 = (1 + k_1 * e ** (-k_2 * (x - x_50))) ** 3
    num_1 = k_1 * k_2 ** 2 * e ** (-k_2 * (x - x_50))
    den_1 = (1 + (k_1 * e ** (-k_2 * (x - x_50)))) ** 2
    v_0 = num_0 / den_0
    v_1 = num_1 / den_1
    return k_3 * (v_0 - v_1)

def calcError(x_0, x_25, x_50, x_75, x_100, k_1, k_2, k_3, k_4):
    """
    Calculate the error of a Cumulative Distribution Function from 
    given x-points and k-transformers.
    Returns the square root of the sum of the distances from CDF
    to the x-points.
    """
    error = sqrt( sum( map( lambda (x,y) : abs(x-y)
        , zip( map( lambda x: CDF(x_0, x_50, x_100, k_1, k_2, k_3, k_4, x)
                , [x_0, x_25, x_50, x_75, x_100] )
                , [0.0, 0.25, 0.50, 0.75, 1.0] ))))
    return error

def calcK2(x_25, x_50, x_75, k_1):
    """
    Calculate k_2 value, given a k_1, and x-points.
    This takes the average of the inner-most x-points, in order
    to make a more smooth CDF.
    """
    if k_1 <=0: return 0
    k_2_1 = -log(k_1/3) / (x_75 - x_50)
    k_2_2 = -log(k_1/3) / (x_50 - x_25)
    k_2   = (k_2_1 + k_2_2) / 2
    return k_2

def calcK3_4(x_0, x_50, x_100, k_1, k_2):
    """
    This adjusts the ends so that they meet at the proper 
    utility value at end points.
    """
    dnom_0 = 1 + (k_1 * e ** (-k_2 * (x_0   - x_50)))
    dnom_1 = 1 + (k_1 * e ** (-k_2 * (x_100 - x_50)))
    
    if dnom_0-dnom_1 == 0: return 0.0, 0.0
    k_3 = (dnom_1 * dnom_0) / (dnom_0 - dnom_1)
    k_4 = (-k_3) / dnom_0

    return k_3, k_4


def GenerateCDFConstructor(x_0, x_25, x_50, x_75, x_100, iterations=13):
    """
    Find the parameters for a CDF, using an iterative method.
    This function finds a fit for the given x-points.
    It returns the  first seven arguments to CDF as a tuple.

    In general the CDF we are building is this:

    CDF(x) = k_3 / (1 + (k_1 * e ** (-k_2 * (x - x_50)))) + k_4


                           k_3
    CDF(x) = -------------------------------- + k_4
                          -k_2 * (x - x_50)
               1 + k_1 * e 

    iterations is a fail-safe, for how many iterations should be done 
    if an steady error cannot be approached.

    My testing shows that this method produces Utility curves that are more
    compacted toward x_50 than they should be.  This algorythm should be 
    replaced with a better fitting algorythm.
    """
    k_1     = 1.0
    k_1_d   = k_1

    k_2 = calcK2(x_25, x_50, x_75, k_1)

    k_1_ret = k_1
    k_2_ret = k_2

    k_3 = 1.0
    k_4 = 0.0

    error = calcError(x_0, x_25, x_50, x_75, x_100, k_1, k_2, k_3, k_4)

    ii, jj = 0, 0
    while ii < iterations:
        k_1_low          = k_1 - k_1_d
        k_2_low          = calcK2(x_25, x_50, x_75, k_1_low)
        k_3_low, k_4_low = calcK3_4(x_0, x_50, x_100, k_1_low, k_2_low)
        error_low        = calcError( x_0, x_25, x_50, x_75, x_100
                                    , k_1_low, k_2_low, k_3_low, k_4_low)

        k_1_high           = k_1 + k_1_d
        k_2_high           = calcK2(x_25, x_50, x_75, k_1_high)
        k_3_high, k_4_high = calcK3_4(x_0, x_50, x_100, k_1_high, k_2_high)
        error_high         = calcError( x_0, x_25, x_50, x_75, x_100 
                                      , k_1_high, k_2_high, k_3_high, k_4_high)


        error_old = error

        if error_high > error_low:
            error = error_low
            k_1   = k_1_low
            k_2   = k_2_low
            k_3   = k_3_low
            k_4   = k_4_low
            ii = 0
        elif error_low > error_high:
            error = error_high
            k_1   = k_1_high
            k_2   = k_2_high
            k_3   = k_3_high
            k_4   = k_4_high
            ii = 0

        if error_old == error : break

        k_1_d = k_1_d / 2
        ii += 1
        jj += 1


    return (x_0, x_50, x_100, k_1, k_2, k_3, k_4)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def combineAngles(xs, ys):
    """
    Create all the combinations of angles from the current lists, and
    the input list.
    xs : [[angles]]
    ys :  [angles]
    """
    result = []
    for x in xs:
        for y in ys:
            newPoints = []
            newPoints.extend(x)
            newPoints.append(y)
            result.append(newPoints)
    return result


def convertHyperSperetoCoordinates(spherePoints):
    """
    spherePoints of the from [r, theta_1, theta_2, theta_3, ... , theta_N]
    a = r * sin(theta_1)
    b = r * cos(theta_1) * sin(theta_2)
    c = r * cos(theta_1) * cos(theta_2) * sin(theta_3)
    d = r * cos(theta_1) * cos(theta_2) * cos(theta_3) * sin(theta_4)
    e = r * cos(theta_1) * cos(theta_2) * cos(theta_3) * cos(theta_4)
    ...
    Note here: 4 angles makes 5 points.
    """
    r      = spherePoints[0]
    theta  = spherePoints[1:]
    points = [r for _ in spherePoints]

    for ii, t in enumerate(theta):
        points[ii] = sin(t) * points[ii]
        points[ii+1:] = map(lambda v: v * cos(t), points[ii+1:])

    return points


def getMSR_ND_hyperShpere_prime( uc_s
                               , r=1.0
                               , N=3
                               , tolerence=1e-9
                               , resultInTheta=False
                               , centerThetas=None
                               , thetaRange=(pi/2.0)
                               ):
    """
    Use a hypsersphere to find the optional progression values
    for the input Utility Curves (uc_s) for a given amount of 
    effort (r) over N evenly spaced test points.
    uc_s : The Utility Curves that are to be optimized.
    r    : The effort (distance from the origin, also radius of the
           hypsersphere) to test the Utility Curves.
    N    : The number of angles to test the hypsersphere at.

    The following arguments are intended for an iterative solution.
    resultInTheta : Returns thetas rather than points.
    centerThetas  : these are the results from less refined iterations
    thetaRange    : the precision of the previous iteration

    Keep in mind the complexity of this function is along O(N^|uc_s|).
    """
    if not centerThetas: centerThetas = [(pi/4.0) for _ in uc_s[1:]]
    thetaStep = thetaRange / float(N)
    
    theta = []
    for centerTheta in centerThetas:
        inner_theta = []
        for ii in range(N+1):
            this_theta = centerTheta + (thetaStep * ((-N/2.0) + ii))
            if   this_theta < 0.0    : this_theta = 0.0
            elif this_theta > pi/2.0 : this_theta = pi/2.0
            inner_theta.append(this_theta)
        theta.append(inner_theta)

    theta[0] = map(lambda x: [r, x], theta[0])

    thetaPairs = reduce(combineAngles, theta)
    points = map(convertHyperSperetoCoordinates, thetaPairs)
    
    values = zip( [ sum( map(lambda (uc, x): uc(x), zip(uc_s, point)) ) 
                    for point in points
                    # if all([x <= uc.x_100 for (uc, x) in zip(uc_s, point)])
                ]
                , points, thetaPairs)

    if len(values) < 1: return None
    maxValue = max(values, key=lambda x: x[0])

    if resultInTheta: return maxValue[2][1:]
    return maxValue[1]


def getMSR_ND_hyperShpere(uc_s, r=1.0, N=3, minStep=1e-12):
    """
    Use a hypsersphere to find the optional progression values
    for the input Utility Curves (uc_s) for a given amount of 
    effort (r) over N evenly spaced test points.
    uc_s    : The Utility Curves that are to be optimized.
    r       : The effort (distance from the origin, also radius of the
              hypsersphere) to test the Utility Curves.
    N       : The number of angles to test the hypsersphere at.
    minStep : The minimum theta range (precision) of the hypsersphere
              sampling.

    Keep in mind the complexity of this function is along O(N^|uc_s|).
    """
    thetas = getMSR_ND_hyperShpere_prime(uc_s, r, N, resultInTheta=True)
    thetaRange  = pi/(2.0 * N)
    while thetaRange > minStep:
        thetas = getMSR_ND_hyperShpere_prime( uc_s, r, N
                                            , resultInTheta=True
                                            , centerThetas=thetas
                                            , thetaRange=thetaRange
                                            )
        thetaRange  = thetaRange / N

    return convertHyperSperetoCoordinates([r] + thetas)
    

def getMSR_ND(uc_s, N=13):
    """
    get the MSR curve for an N-dimensional option set, represented by
    the Utility Curves (uc_s).  The result will be the optimized progression
    for the Utility Curves for a given effort, not the sequence that should
    be followed.

    Keep in mind the complexity of this function is along O(N^(|uc_s|+1)).
    """
    radiusMin = 0.0
    radiusMax = sqrt( sum([uc.denormalize(uc.x_100) ** 2 for uc in uc_s]) )
    radiusStep = (radiusMax - radiusMin) / float(N)

    radii = ( radiusStep * ii for ii in range(N+1) )
    
    points = [ getMSR_ND_hyperShpere(uc_s,r) 
               for r in radii
             ]

    points = filter(lambda x: x!=None, points)
    return points




if __name__=="__main__":
    from pprint import pprint
    A = UtilityCurve(0,5,7,12,20)
    B = UtilityCurve(0,3,4,7,10)
    C = UtilityCurve(10,15,71,112,210)
    D = UtilityCurve(1000,8903,10006,100190,10020)
    E = UtilityCurve(1, 6, 15, 20, 30)

    uc_s = [A,B,C,D,E]
    points = getMSR_ND(uc_s)
    pprint(points)



