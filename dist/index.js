(function (factory) {
    typeof define === 'function' && define.amd ? define(factory) :
    factory();
}(function () { 'use strict';

    function tinyNDArrayOfInteger (gridShape) {
        var dimensions = gridShape.length,
            totalLength = 1,
            stride = new Array(dimensions),
            dimension;

        for (dimension = dimensions; dimension > 0; dimension--) {
            stride[dimension - 1] = totalLength;
            totalLength = totalLength * gridShape[dimension - 1];
        }

        return {
            stride: stride,
            data: new Uint32Array(totalLength)
        };
    }

    function tinyNDArrayOfArray (gridShape) {
        var dimensions = gridShape.length,
            totalLength = 1,
            stride = new Array(dimensions),
            data = [],
            dimension, index;

        for (dimension = dimensions; dimension > 0; dimension--) {
            stride[dimension - 1] = totalLength;
            totalLength = totalLength * gridShape[dimension - 1];
        }

        for (index = 0; index < totalLength; index++) {
            data.push([]);
        }

        return {
            stride: stride,
            data: data
        };
    }

    var tinyNdarray = {
        integer: tinyNDArrayOfInteger,
        array: tinyNDArrayOfArray
    };

    // sphere-random module by Mikola Lysenko under the MIT License
    // waiting for https://github.com/scijs/sphere-random/pull/1 to be merged

    var sphereRandom = sampleSphere;

    /**
     * @param {int} d Dimensions
     * @param {Function} rng
     * @returns {Array}
     */
    function sampleSphere(d, rng) {
        var v = new Array(d),
            d2 = Math.floor(d/2) << 1,
            r2 = 0.0,
            rr,
            r,
            theta,
            h,
            i;

        for (i = 0; i < d2; i += 2) {
            rr = -2.0 * Math.log(rng());
            r =  Math.sqrt(rr);
            theta = 2.0 * Math.PI * rng();

            r2+= rr;
            v[i] = r * Math.cos(theta);
            v[i+1] = r * Math.sin(theta);
        }

        if (d % 2) {
            var x = Math.sqrt(-2.0 * Math.log(rng())) * Math.cos(2.0 * Math.PI * rng());
            v[d - 1] = x;
            r2+= Math.pow(x, 2);
        }

        h = 1.0 / Math.sqrt(r2);

        for (i = 0; i < d; ++i) {
            v[i] *= h;
        }

        return v;
    }

    var moore = function moore(range, dimensions) {
      range = range || 1;
      dimensions = dimensions || 2;

      var size = range * 2 + 1;
      var length = Math.pow(size, dimensions) - 1;
      var neighbors = new Array(length);

      for (var i = 0; i < length; i++) {
        var neighbor = neighbors[i] = new Array(dimensions);
        var index = i < length / 2 ? i : i + 1;
        for (var dimension = 1; dimension <= dimensions; dimension++) {
          var value = index % Math.pow(size, dimension);
          neighbor[dimension - 1] = value / Math.pow(size, dimension - 1) - range;
          index -= value;
        }
      }

      return neighbors
    };

    /**
     * Get the neighbourhood ordered by distance, including the origin point
     * @param {int} dimensionNumber Number of dimensions
     * @returns {Array} Neighbourhood
     */
    function getNeighbourhood (dimensionNumber) {
        var neighbourhood = moore(2, dimensionNumber),
            origin = [],
            dimension;

        for (dimension = 0; dimension < dimensionNumber; dimension++) {
            origin.push(0);
        }

        neighbourhood.push(origin);

        // sort by ascending distance to optimize proximity checks
        // see point 5.1 in Parallel Poisson Disk Sampling by Li-Yi Wei, 2008
        // http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.460.3061&rank=1
        neighbourhood.sort(function (n1, n2) {
            var squareDist1 = 0,
                squareDist2 = 0,
                dimension;

            for (dimension = 0; dimension < dimensionNumber; dimension++) {
                squareDist1 += Math.pow(n1[dimension], 2);
                squareDist2 += Math.pow(n2[dimension], 2);
            }

            if (squareDist1 < squareDist2) {
                return -1;
            } else if(squareDist1 > squareDist2) {
                return 1;
            } else {
                return 0;
            }
        });

        return neighbourhood;
    }

    var neighbourhood = getNeighbourhood;

    var tinyNDArray = tinyNdarray.integer;

    /**
     * Get the squared euclidean distance from two points of arbitrary, but equal, dimensions
     * @param {Array} point1
     * @param {Array} point2
     * @returns {number} Squared euclidean distance
     */
    function squaredEuclideanDistance (point1, point2) {
        var result = 0,
            i = 0;

        for (; i < point1.length; i++) {
            result += Math.pow(point1[i] - point2[i], 2);
        }

        return result;
    }

    /**
     * FixedDensityPDS constructor
     * @param {object} options Options
     * @param {Array} options.shape Shape of the space
     * @param {float} options.minDistance Minimum distance between each points
     * @param {float} [options.maxDistance] Maximum distance between each points, defaults to minDistance * 2
     * @param {int} [options.tries] Number of times the algorithm will try to place a point in the neighbourhood of another points, defaults to 30
     * @param {function|null} [rng] RNG function, defaults to Math.random
     * @constructor
     */
    function FixedDensityPDS (options, rng) {
        if (typeof options.distanceFunction === 'function') {
            throw new Error('PoissonDiskSampling: Tried to instantiate the fixed density implementation with a distanceFunction');
        }

        this.shape = options.shape;
        this.minDistance = options.minDistance;
        this.maxDistance = options.maxDistance || options.minDistance * 2;
        this.maxTries = Math.ceil(Math.max(1, options.tries || 30));

        this.rng = rng || Math.random;

        this.dimension = this.shape.length;
        this.squaredMinDistance = this.minDistance * this.minDistance;
        this.deltaDistance = this.maxDistance - this.minDistance;
        this.cellSize = this.minDistance / Math.sqrt(this.dimension);

        this.neighbourhood = neighbourhood(this.dimension);

        this.currentPoint = null;
        this.processList = [];
        this.samplePoints = [];

        // cache grid

        this.gridShape = [];

        for (var i = 0; i < this.dimension; i++) {
            this.gridShape.push(Math.ceil(this.shape[i] / this.cellSize));
        }

        this.grid = tinyNDArray(this.gridShape); //will store references to samplePoints
    }

    FixedDensityPDS.prototype.shape = null;
    FixedDensityPDS.prototype.dimension = null;
    FixedDensityPDS.prototype.minDistance = null;
    FixedDensityPDS.prototype.maxDistance = null;
    FixedDensityPDS.prototype.squaredMinDistance = null;
    FixedDensityPDS.prototype.deltaDistance = null;
    FixedDensityPDS.prototype.cellSize = null;
    FixedDensityPDS.prototype.maxTries = null;
    FixedDensityPDS.prototype.rng = null;
    FixedDensityPDS.prototype.neighbourhood = null;

    FixedDensityPDS.prototype.currentPoint = null;
    FixedDensityPDS.prototype.processList = null;
    FixedDensityPDS.prototype.samplePoints = null;
    FixedDensityPDS.prototype.gridShape = null;
    FixedDensityPDS.prototype.grid = null;

    /**
     * Add a totally random point in the grid
     * @returns {Array} The point added to the grid
     */
    FixedDensityPDS.prototype.addRandomPoint = function () {
        var point = new Array(this.dimension);

        for (var i = 0; i < this.dimension; i++) {
            point[i] = this.rng() * this.shape[i];
        }

        return this.directAddPoint(point);
    };

    /**
     * Add a given point to the grid
     * @param {Array} point Point
     * @returns {Array|null} The point added to the grid, null if the point is out of the bound or not of the correct dimension
     */
    FixedDensityPDS.prototype.addPoint = function (point) {
        var dimension,
            valid = true;

        if (point.length === this.dimension) {
            for (dimension = 0; dimension < this.dimension && valid; dimension++) {
                valid = (point[dimension] >= 0 && point[dimension] <= this.shape[dimension]);
            }
        } else {
            valid = false;
        }

        return valid ? this.directAddPoint(point) : null;
    };

    /**
     * Add a given point to the grid, without any check
     * @param {Array} point Point
     * @returns {Array} The point added to the grid
     * @protected
     */
    FixedDensityPDS.prototype.directAddPoint = function (point) {
        var internalArrayIndex = 0,
            stride = this.grid.stride,
            dimension;

        this.processList.push(point);
        this.samplePoints.push(point);

        for (dimension = 0; dimension < this.dimension; dimension++) {
            internalArrayIndex += ((point[dimension] / this.cellSize) | 0) * stride[dimension];
        }

        this.grid.data[internalArrayIndex] = this.samplePoints.length; // store the point reference

        return point;
    };

    /**
     * Check whether a given point is in the neighbourhood of existing points
     * @param {Array} point Point
     * @returns {boolean} Whether the point is in the neighbourhood of another point
     * @protected
     */
    FixedDensityPDS.prototype.inNeighbourhood = function (point) {
        var dimensionNumber = this.dimension,
            stride = this.grid.stride,
            neighbourIndex,
            internalArrayIndex,
            dimension,
            currentDimensionValue,
            existingPoint;

        for (neighbourIndex = 0; neighbourIndex < this.neighbourhood.length; neighbourIndex++) {
            internalArrayIndex = 0;

            for (dimension = 0; dimension < dimensionNumber; dimension++) {
                currentDimensionValue = ((point[dimension] / this.cellSize) | 0) + this.neighbourhood[neighbourIndex][dimension];

                if (currentDimensionValue >= 0 && currentDimensionValue < this.gridShape[dimension]) {
                    internalArrayIndex += currentDimensionValue * stride[dimension];
                }
            }

            if (this.grid.data[internalArrayIndex] !== 0) {
                existingPoint = this.samplePoints[this.grid.data[internalArrayIndex] - 1];

                if (squaredEuclideanDistance(point, existingPoint) < this.squaredMinDistance) {
                    return true;
                }
            }
        }

        return false;
    };

    /**
     * Try to generate a new point in the grid, returns null if it wasn't possible
     * @returns {Array|null} The added point or null
     */
    FixedDensityPDS.prototype.next = function () {
        var tries,
            angle,
            distance,
            currentPoint,
            newPoint,
            inShape,
            i;

        while (this.processList.length > 0) {
            if (this.currentPoint === null) {
                this.currentPoint = this.processList.shift();
            }

            currentPoint = this.currentPoint;

            for (tries = 0; tries < this.maxTries; tries++) {
                inShape = true;
                distance = this.minDistance + this.deltaDistance * this.rng();

                if (this.dimension === 2) {
                    angle = this.rng() * Math.PI * 2;
                    newPoint = [
                        Math.cos(angle),
                        Math.sin(angle)
                    ];
                } else {
                    newPoint = sphereRandom(this.dimension, this.rng);
                }

                for (i = 0; inShape && i < this.dimension; i++) {
                    newPoint[i] = currentPoint[i] + newPoint[i] * distance;
                    inShape = (newPoint[i] >= 0 && newPoint[i] <= this.shape[i] - 1);
                }

                if (inShape && !this.inNeighbourhood(newPoint)) {
                    return this.directAddPoint(newPoint);
                }
            }

            if (tries === this.maxTries) {
                this.currentPoint = null;
            }
        }

        return null;
    };

    /**
     * Automatically fill the grid, adding a random point to start the process if needed.
     * Will block the thread, probably best to use it in a web worker or child process.
     * @returns {Array[]} Sample points
     */
    FixedDensityPDS.prototype.fill = function () {
        if (this.samplePoints.length === 0) {
            this.addRandomPoint();
        }

        while(this.next()) {}

        return this.samplePoints;
    };

    /**
     * Get all the points in the grid.
     * @returns {Array[]} Sample points
     */
    FixedDensityPDS.prototype.getAllPoints = function () {
        return this.samplePoints;
    };

    /**
     * Get all the points in the grid along with the result of the distance function.
     * @throws Will always throw an error.
     */
    FixedDensityPDS.prototype.getAllPointsWithDistance = function () {
        throw new Error('PoissonDiskSampling: getAllPointsWithDistance() is not available in fixed-density implementation');
    };

    /**
     * Reinitialize the grid as well as the internal state
     */
    FixedDensityPDS.prototype.reset = function () {
        var gridData = this.grid.data,
            i = 0;

        // reset the cache grid
        for (i = 0; i < gridData.length; i++) {
            gridData[i] = 0;
        }

        // new array for the samplePoints as it is passed by reference to the outside
        this.samplePoints = [];

        // reset the internal state
        this.currentPoint = null;
        this.processList.length = 0;
    };

    var fixedDensity = FixedDensityPDS;

    var tinyNDArray$1 = tinyNdarray.array;

    /**
     * Get the euclidean distance from two points of arbitrary, but equal, dimensions
     * @param {Array} point1
     * @param {Array} point2
     * @returns {number} Euclidean distance
     */
    function euclideanDistance (point1, point2) {
        var result = 0,
            i = 0;

        for (; i < point1.length; i++) {
            result += Math.pow(point1[i] - point2[i], 2);
        }

        return Math.sqrt(result);
    }

    /**
     * VariableDensityPDS constructor
     * @param {object} options Options
     * @param {Array} options.shape Shape of the space
     * @param {float} options.minDistance Minimum distance between each points
     * @param {float} [options.maxDistance] Maximum distance between each points, defaults to minDistance * 2
     * @param {int} [options.tries] Number of times the algorithm will try to place a point in the neighbourhood of another points, defaults to 30
     * @param {function} options.distanceFunction Function to control the distance between each point depending on their position, must return a value between 0 and 1
     * @param {float} [options.bias] When using a distanceFunction, will indicate which point constraint takes priority when evaluating two points (0 for the lowest distance, 1 for the highest distance), defaults to 0
     * @param {function|null} rng RNG function, defaults to Math.random
     * @constructor
     */
    function VariableDensityPDS (options, rng) {
        if (typeof options.distanceFunction !== 'function') {
            throw new Error('PoissonDiskSampling: Tried to instantiate the variable density implementation without a distanceFunction');
        }

        this.shape = options.shape;
        this.minDistance = options.minDistance;
        this.maxDistance = options.maxDistance || options.minDistance * 2;
        this.maxTries = Math.ceil(Math.max(1, options.tries || 30));
        this.distanceFunction = options.distanceFunction;
        this.bias = Math.max(0, Math.min(1, options.bias || 0));

        this.rng = rng || Math.random;

        this.dimension = this.shape.length;
        this.deltaDistance = this.maxDistance - this.minDistance;
        this.cellSize = this.maxDistance / Math.sqrt(this.dimension);

        this.neighbourhood = neighbourhood(this.dimension);

        this.currentPoint = null;
        this.currentDistance = 0;
        this.processList = [];
        this.samplePoints = [];
        this.sampleDistance = []; // used to store the distance for a given point

        // cache grid

        this.gridShape = [];

        for (var i = 0; i < this.dimension; i++) {
            this.gridShape.push(Math.ceil(this.shape[i] / this.cellSize));
        }

        this.grid = tinyNDArray$1(this.gridShape); //will store references to samplePoints and sampleDistance
    }

    VariableDensityPDS.prototype.shape = null;
    VariableDensityPDS.prototype.dimension = null;
    VariableDensityPDS.prototype.minDistance = null;
    VariableDensityPDS.prototype.maxDistance = null;
    VariableDensityPDS.prototype.deltaDistance = null;
    VariableDensityPDS.prototype.cellSize = null;
    VariableDensityPDS.prototype.maxTries = null;
    VariableDensityPDS.prototype.distanceFunction = null;
    VariableDensityPDS.prototype.bias = null;
    VariableDensityPDS.prototype.rng = null;
    VariableDensityPDS.prototype.neighbourhood = null;

    VariableDensityPDS.prototype.currentPoint = null;
    VariableDensityPDS.prototype.currentDistance = null;
    VariableDensityPDS.prototype.processList = null;
    VariableDensityPDS.prototype.samplePoints = null;
    VariableDensityPDS.prototype.sampleDistance = null;
    VariableDensityPDS.prototype.gridShape = null;
    VariableDensityPDS.prototype.grid = null;

    /**
     * Add a totally random point in the grid
     * @returns {Array} The point added to the grid
     */
    VariableDensityPDS.prototype.addRandomPoint = function () {
        var point = new Array(this.dimension);

        for (var i = 0; i < this.dimension; i++) {
            point[i] = this.rng() * this.shape[i];
        }

        return this.directAddPoint(point);
    };

    /**
     * Add a given point to the grid
     * @param {Array} point Point
     * @returns {Array|null} The point added to the grid, null if the point is out of the bound or not of the correct dimension
     */
    VariableDensityPDS.prototype.addPoint = function (point) {
        var dimension,
            valid = true;

        if (point.length === this.dimension) {
            for (dimension = 0; dimension < this.dimension && valid; dimension++) {
                valid = (point[dimension] >= 0 && point[dimension] <= this.shape[dimension]);
            }
        } else {
            valid = false;
        }

        return valid ? this.directAddPoint(point) : null;
    };

    /**
     * Add a given point to the grid, without any check
     * @param {Array} point Point
     * @returns {Array} The point added to the grid
     * @protected
     */
    VariableDensityPDS.prototype.directAddPoint = function (point) {
        var internalArrayIndex = 0,
            stride = this.grid.stride,
            pointIndex = this.samplePoints.length,
            dimension;

        this.processList.push(pointIndex);
        this.samplePoints.push(point);
        this.sampleDistance.push(this.distanceFunction(point));

        for (dimension = 0; dimension < this.dimension; dimension++) {
            internalArrayIndex += ((point[dimension] / this.cellSize) | 0) * stride[dimension];
        }

        this.grid.data[internalArrayIndex].push(pointIndex); // store the point reference

        return point;
    };

    /**
     * Check whether a given point is in the neighbourhood of existing points
     * @param {Array} point Point
     * @returns {boolean} Whether the point is in the neighbourhood of another point
     * @protected
     */
    VariableDensityPDS.prototype.inNeighbourhood = function (point) {
        var dimensionNumber = this.dimension,
            stride = this.grid.stride,
            neighbourIndex,
            internalArrayIndex,
            dimension,
            currentDimensionValue,
            existingPoint,
            existingPointDistance;

        var pointDistance = this.distanceFunction(point);

        for (neighbourIndex = 0; neighbourIndex < this.neighbourhood.length; neighbourIndex++) {
            internalArrayIndex = 0;

            for (dimension = 0; dimension < dimensionNumber; dimension++) {
                currentDimensionValue = ((point[dimension] / this.cellSize) | 0) + this.neighbourhood[neighbourIndex][dimension];

                if (currentDimensionValue >= 0 && currentDimensionValue < this.gridShape[dimension]) {
                    internalArrayIndex += currentDimensionValue * stride[dimension];
                }
            }

            if (this.grid.data[internalArrayIndex].length > 0) {
                for (var i = 0; i < this.grid.data[internalArrayIndex].length; i++) {
                    existingPoint = this.samplePoints[this.grid.data[internalArrayIndex][i]];
                    existingPointDistance = this.sampleDistance[this.grid.data[internalArrayIndex][i]];

                    var minDistance = Math.min(existingPointDistance, pointDistance);
                    var maxDistance = Math.max(existingPointDistance, pointDistance);
                    var dist = minDistance + (maxDistance - minDistance) * this.bias;


                    if (euclideanDistance(point, existingPoint) < this.minDistance + this.deltaDistance * dist) {
                        return true;
                    }
                }
            }
        }

        return false;
    };

    /**
     * Try to generate a new point in the grid, returns null if it wasn't possible
     * @returns {Array|null} The added point or null
     */
    VariableDensityPDS.prototype.next = function () {
        var tries,
            angle,
            distance,
            currentPoint,
            currentDistance,
            newPoint,
            inShape,
            i;

        while (this.processList.length > 0) {
            if (this.currentPoint === null) {
                var sampleIndex = this.processList.shift();
                this.currentPoint = this.samplePoints[sampleIndex];
                this.currentDistance = this.sampleDistance[sampleIndex];
            }

            currentPoint = this.currentPoint;
            currentDistance = this.currentDistance;

            for (tries = 0; tries < this.maxTries; tries++) {
                inShape = true;
                distance = this.minDistance + this.deltaDistance * (currentDistance + (1 - currentDistance) * this.bias);

                if (this.dimension === 2) {
                    angle = this.rng() * Math.PI * 2;
                    newPoint = [
                        Math.cos(angle),
                        Math.sin(angle)
                    ];
                } else {
                    newPoint = sphereRandom(this.dimension, this.rng);
                }

                for (i = 0; inShape && i < this.dimension; i++) {
                    newPoint[i] = currentPoint[i] + newPoint[i] * distance;
                    inShape = (newPoint[i] >= 0 && newPoint[i] <= this.shape[i] - 1);
                }

                if (inShape && !this.inNeighbourhood(newPoint)) {
                    return this.directAddPoint(newPoint);
                }
            }

            if (tries === this.maxTries) {
                this.currentPoint = null;
            }
        }

        return null;
    };

    /**
     * Automatically fill the grid, adding a random point to start the process if needed.
     * Will block the thread, probably best to use it in a web worker or child process.
     * @returns {Array[]} Sample points
     */
    VariableDensityPDS.prototype.fill = function () {
        if (this.samplePoints.length === 0) {
            this.addRandomPoint();
        }

        while(this.next()) {}

        return this.samplePoints;
    };

    /**
     * Get all the points in the grid.
     * @returns {Array[]} Sample points
     */
    VariableDensityPDS.prototype.getAllPoints = function () {
        return this.samplePoints;
    };

    /**
     * Get all the points in the grid along with the result of the distance function.
     * @returns {Array[]} Sample points with their distance function result
     */
    VariableDensityPDS.prototype.getAllPointsWithDistance = function () {
        var result = new Array(this.samplePoints.length),
            i = 0,
            dimension = 0,
            point;

        for (i = 0; i < this.samplePoints.length; i++) {
            point = new Array(this.dimension + 1);

            for (dimension = 0; dimension < this.dimension; dimension++) {
                point[dimension] = this.samplePoints[i][dimension];
            }

            point[this.dimension] = this.sampleDistance[i];

            result[i] = point;
        }

        return result;
    };

    /**
     * Reinitialize the grid as well as the internal state
     */
    VariableDensityPDS.prototype.reset = function () {
        var gridData = this.grid.data,
            i = 0;

        // reset the cache grid
        for (i = 0; i < gridData.length; i++) {
            gridData[i] = [];
        }

        // new array for the samplePoints as it is passed by reference to the outside
        this.samplePoints = [];

        // reset the internal state
        this.currentPoint = null;
        this.processList.length = 0;
    };

    var variableDensity = VariableDensityPDS;

    /**
     * PoissonDiskSampling constructor
     * @param {object} options Options
     * @param {Array} options.shape Shape of the space
     * @param {float} options.minDistance Minimum distance between each points
     * @param {float} [options.maxDistance] Maximum distance between each points, defaults to minDistance * 2
     * @param {int} [options.tries] Number of times the algorithm will try to place a point in the neighbourhood of another points, defaults to 30
     * @param {function|null} [options.distanceFunction] Function to control the distance between each point depending on their position, must return a value between 0 and 1
     * @param {function|null} [options.bias] When using a distanceFunction, will indicate which point constraint takes priority when evaluating two points (0 for the lowest distance, 1 for the highest distance), defaults to 0
     * @param {function|null} [rng] RNG function, defaults to Math.random
     * @constructor
     */
    function PoissonDiskSampling (options, rng) {
        this.shape = options.shape;

        if (typeof options.distanceFunction === 'function') {
            this.implementation = new variableDensity(options, rng);
        } else {
            this.implementation = new fixedDensity(options, rng);
        }
    }

    PoissonDiskSampling.prototype.implementation = null;

    /**
     * Add a totally random point in the grid
     * @returns {Array} The point added to the grid
     */
    PoissonDiskSampling.prototype.addRandomPoint = function () {
        return this.implementation.addRandomPoint();
    };

    /**
     * Add a given point to the grid
     * @param {Array} point Point
     * @returns {Array|null} The point added to the grid, null if the point is out of the bound or not of the correct dimension
     */
    PoissonDiskSampling.prototype.addPoint = function (point) {
        return this.implementation.addPoint(point);
    };

    /**
     * Try to generate a new point in the grid, returns null if it wasn't possible
     * @returns {Array|null} The added point or null
     */
    PoissonDiskSampling.prototype.next = function () {
        return this.implementation.next();
    };

    /**
     * Automatically fill the grid, adding a random point to start the process if needed.
     * Will block the thread, probably best to use it in a web worker or child process.
     * @returns {Array[]} Sample points
     */
    PoissonDiskSampling.prototype.fill = function () {
        return this.implementation.fill();
    };

    /**
     * Get all the points in the grid.
     * @returns {Array[]} Sample points
     */
    PoissonDiskSampling.prototype.getAllPoints = function () {
        return this.implementation.getAllPoints();
    };

    /**
     * Get all the points in the grid along with the result of the distance function.
     * @throws Will throw an error if a distance function was not provided to the constructor.
     * @returns {Array[]} Sample points with their distance function result
     */
    PoissonDiskSampling.prototype.getAllPointsWithDistance = function () {
        return this.implementation.getAllPointsWithDistance();
    };

    /**
     * Reinitialize the grid as well as the internal state
     */
    PoissonDiskSampling.prototype.reset = function () {
        this.implementation.reset();
    };

    var poissonDiskSampling = PoissonDiskSampling;

    function createCommonjsModule(fn, module) {
    	return module = { exports: {} }, fn(module, module.exports), module.exports;
    }

    var simplexNoise = createCommonjsModule(function (module, exports) {
    /*
     * A fast javascript implementation of simplex noise by Jonas Wagner

    Based on a speed-improved simplex noise algorithm for 2D, 3D and 4D in Java.
    Which is based on example code by Stefan Gustavson (stegu@itn.liu.se).
    With Optimisations by Peter Eastman (peastman@drizzle.stanford.edu).
    Better rank ordering method by Stefan Gustavson in 2012.


     Copyright (c) 2018 Jonas Wagner

     Permission is hereby granted, free of charge, to any person obtaining a copy
     of this software and associated documentation files (the "Software"), to deal
     in the Software without restriction, including without limitation the rights
     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
     copies of the Software, and to permit persons to whom the Software is
     furnished to do so, subject to the following conditions:

     The above copyright notice and this permission notice shall be included in all
     copies or substantial portions of the Software.

     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
     SOFTWARE.
     */
    (function() {

      var F2 = 0.5 * (Math.sqrt(3.0) - 1.0);
      var G2 = (3.0 - Math.sqrt(3.0)) / 6.0;
      var F3 = 1.0 / 3.0;
      var G3 = 1.0 / 6.0;
      var F4 = (Math.sqrt(5.0) - 1.0) / 4.0;
      var G4 = (5.0 - Math.sqrt(5.0)) / 20.0;

      function SimplexNoise(randomOrSeed) {
        var random;
        if (typeof randomOrSeed == 'function') {
          random = randomOrSeed;
        }
        else if (randomOrSeed) {
          random = alea(randomOrSeed);
        } else {
          random = Math.random;
        }
        this.p = buildPermutationTable(random);
        this.perm = new Uint8Array(512);
        this.permMod12 = new Uint8Array(512);
        for (var i = 0; i < 512; i++) {
          this.perm[i] = this.p[i & 255];
          this.permMod12[i] = this.perm[i] % 12;
        }

      }
      SimplexNoise.prototype = {
        grad3: new Float32Array([1, 1, 0,
          -1, 1, 0,
          1, -1, 0,

          -1, -1, 0,
          1, 0, 1,
          -1, 0, 1,

          1, 0, -1,
          -1, 0, -1,
          0, 1, 1,

          0, -1, 1,
          0, 1, -1,
          0, -1, -1]),
        grad4: new Float32Array([0, 1, 1, 1, 0, 1, 1, -1, 0, 1, -1, 1, 0, 1, -1, -1,
          0, -1, 1, 1, 0, -1, 1, -1, 0, -1, -1, 1, 0, -1, -1, -1,
          1, 0, 1, 1, 1, 0, 1, -1, 1, 0, -1, 1, 1, 0, -1, -1,
          -1, 0, 1, 1, -1, 0, 1, -1, -1, 0, -1, 1, -1, 0, -1, -1,
          1, 1, 0, 1, 1, 1, 0, -1, 1, -1, 0, 1, 1, -1, 0, -1,
          -1, 1, 0, 1, -1, 1, 0, -1, -1, -1, 0, 1, -1, -1, 0, -1,
          1, 1, 1, 0, 1, 1, -1, 0, 1, -1, 1, 0, 1, -1, -1, 0,
          -1, 1, 1, 0, -1, 1, -1, 0, -1, -1, 1, 0, -1, -1, -1, 0]),
        noise2D: function(xin, yin) {
          var permMod12 = this.permMod12;
          var perm = this.perm;
          var grad3 = this.grad3;
          var n0 = 0; // Noise contributions from the three corners
          var n1 = 0;
          var n2 = 0;
          // Skew the input space to determine which simplex cell we're in
          var s = (xin + yin) * F2; // Hairy factor for 2D
          var i = Math.floor(xin + s);
          var j = Math.floor(yin + s);
          var t = (i + j) * G2;
          var X0 = i - t; // Unskew the cell origin back to (x,y) space
          var Y0 = j - t;
          var x0 = xin - X0; // The x,y distances from the cell origin
          var y0 = yin - Y0;
          // For the 2D case, the simplex shape is an equilateral triangle.
          // Determine which simplex we are in.
          var i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
          if (x0 > y0) {
            i1 = 1;
            j1 = 0;
          } // lower triangle, XY order: (0,0)->(1,0)->(1,1)
          else {
            i1 = 0;
            j1 = 1;
          } // upper triangle, YX order: (0,0)->(0,1)->(1,1)
          // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
          // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
          // c = (3-sqrt(3))/6
          var x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
          var y1 = y0 - j1 + G2;
          var x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
          var y2 = y0 - 1.0 + 2.0 * G2;
          // Work out the hashed gradient indices of the three simplex corners
          var ii = i & 255;
          var jj = j & 255;
          // Calculate the contribution from the three corners
          var t0 = 0.5 - x0 * x0 - y0 * y0;
          if (t0 >= 0) {
            var gi0 = permMod12[ii + perm[jj]] * 3;
            t0 *= t0;
            n0 = t0 * t0 * (grad3[gi0] * x0 + grad3[gi0 + 1] * y0); // (x,y) of grad3 used for 2D gradient
          }
          var t1 = 0.5 - x1 * x1 - y1 * y1;
          if (t1 >= 0) {
            var gi1 = permMod12[ii + i1 + perm[jj + j1]] * 3;
            t1 *= t1;
            n1 = t1 * t1 * (grad3[gi1] * x1 + grad3[gi1 + 1] * y1);
          }
          var t2 = 0.5 - x2 * x2 - y2 * y2;
          if (t2 >= 0) {
            var gi2 = permMod12[ii + 1 + perm[jj + 1]] * 3;
            t2 *= t2;
            n2 = t2 * t2 * (grad3[gi2] * x2 + grad3[gi2 + 1] * y2);
          }
          // Add contributions from each corner to get the final noise value.
          // The result is scaled to return values in the interval [-1,1].
          return 70.0 * (n0 + n1 + n2);
        },
        // 3D simplex noise
        noise3D: function(xin, yin, zin) {
          var permMod12 = this.permMod12;
          var perm = this.perm;
          var grad3 = this.grad3;
          var n0, n1, n2, n3; // Noise contributions from the four corners
          // Skew the input space to determine which simplex cell we're in
          var s = (xin + yin + zin) * F3; // Very nice and simple skew factor for 3D
          var i = Math.floor(xin + s);
          var j = Math.floor(yin + s);
          var k = Math.floor(zin + s);
          var t = (i + j + k) * G3;
          var X0 = i - t; // Unskew the cell origin back to (x,y,z) space
          var Y0 = j - t;
          var Z0 = k - t;
          var x0 = xin - X0; // The x,y,z distances from the cell origin
          var y0 = yin - Y0;
          var z0 = zin - Z0;
          // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
          // Determine which simplex we are in.
          var i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
          var i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
          if (x0 >= y0) {
            if (y0 >= z0) {
              i1 = 1;
              j1 = 0;
              k1 = 0;
              i2 = 1;
              j2 = 1;
              k2 = 0;
            } // X Y Z order
            else if (x0 >= z0) {
              i1 = 1;
              j1 = 0;
              k1 = 0;
              i2 = 1;
              j2 = 0;
              k2 = 1;
            } // X Z Y order
            else {
              i1 = 0;
              j1 = 0;
              k1 = 1;
              i2 = 1;
              j2 = 0;
              k2 = 1;
            } // Z X Y order
          }
          else { // x0<y0
            if (y0 < z0) {
              i1 = 0;
              j1 = 0;
              k1 = 1;
              i2 = 0;
              j2 = 1;
              k2 = 1;
            } // Z Y X order
            else if (x0 < z0) {
              i1 = 0;
              j1 = 1;
              k1 = 0;
              i2 = 0;
              j2 = 1;
              k2 = 1;
            } // Y Z X order
            else {
              i1 = 0;
              j1 = 1;
              k1 = 0;
              i2 = 1;
              j2 = 1;
              k2 = 0;
            } // Y X Z order
          }
          // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
          // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
          // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
          // c = 1/6.
          var x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
          var y1 = y0 - j1 + G3;
          var z1 = z0 - k1 + G3;
          var x2 = x0 - i2 + 2.0 * G3; // Offsets for third corner in (x,y,z) coords
          var y2 = y0 - j2 + 2.0 * G3;
          var z2 = z0 - k2 + 2.0 * G3;
          var x3 = x0 - 1.0 + 3.0 * G3; // Offsets for last corner in (x,y,z) coords
          var y3 = y0 - 1.0 + 3.0 * G3;
          var z3 = z0 - 1.0 + 3.0 * G3;
          // Work out the hashed gradient indices of the four simplex corners
          var ii = i & 255;
          var jj = j & 255;
          var kk = k & 255;
          // Calculate the contribution from the four corners
          var t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
          if (t0 < 0) n0 = 0.0;
          else {
            var gi0 = permMod12[ii + perm[jj + perm[kk]]] * 3;
            t0 *= t0;
            n0 = t0 * t0 * (grad3[gi0] * x0 + grad3[gi0 + 1] * y0 + grad3[gi0 + 2] * z0);
          }
          var t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
          if (t1 < 0) n1 = 0.0;
          else {
            var gi1 = permMod12[ii + i1 + perm[jj + j1 + perm[kk + k1]]] * 3;
            t1 *= t1;
            n1 = t1 * t1 * (grad3[gi1] * x1 + grad3[gi1 + 1] * y1 + grad3[gi1 + 2] * z1);
          }
          var t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
          if (t2 < 0) n2 = 0.0;
          else {
            var gi2 = permMod12[ii + i2 + perm[jj + j2 + perm[kk + k2]]] * 3;
            t2 *= t2;
            n2 = t2 * t2 * (grad3[gi2] * x2 + grad3[gi2 + 1] * y2 + grad3[gi2 + 2] * z2);
          }
          var t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
          if (t3 < 0) n3 = 0.0;
          else {
            var gi3 = permMod12[ii + 1 + perm[jj + 1 + perm[kk + 1]]] * 3;
            t3 *= t3;
            n3 = t3 * t3 * (grad3[gi3] * x3 + grad3[gi3 + 1] * y3 + grad3[gi3 + 2] * z3);
          }
          // Add contributions from each corner to get the final noise value.
          // The result is scaled to stay just inside [-1,1]
          return 32.0 * (n0 + n1 + n2 + n3);
        },
        // 4D simplex noise, better simplex rank ordering method 2012-03-09
        noise4D: function(x, y, z, w) {
          var perm = this.perm;
          var grad4 = this.grad4;

          var n0, n1, n2, n3, n4; // Noise contributions from the five corners
          // Skew the (x,y,z,w) space to determine which cell of 24 simplices we're in
          var s = (x + y + z + w) * F4; // Factor for 4D skewing
          var i = Math.floor(x + s);
          var j = Math.floor(y + s);
          var k = Math.floor(z + s);
          var l = Math.floor(w + s);
          var t = (i + j + k + l) * G4; // Factor for 4D unskewing
          var X0 = i - t; // Unskew the cell origin back to (x,y,z,w) space
          var Y0 = j - t;
          var Z0 = k - t;
          var W0 = l - t;
          var x0 = x - X0; // The x,y,z,w distances from the cell origin
          var y0 = y - Y0;
          var z0 = z - Z0;
          var w0 = w - W0;
          // For the 4D case, the simplex is a 4D shape I won't even try to describe.
          // To find out which of the 24 possible simplices we're in, we need to
          // determine the magnitude ordering of x0, y0, z0 and w0.
          // Six pair-wise comparisons are performed between each possible pair
          // of the four coordinates, and the results are used to rank the numbers.
          var rankx = 0;
          var ranky = 0;
          var rankz = 0;
          var rankw = 0;
          if (x0 > y0) rankx++;
          else ranky++;
          if (x0 > z0) rankx++;
          else rankz++;
          if (x0 > w0) rankx++;
          else rankw++;
          if (y0 > z0) ranky++;
          else rankz++;
          if (y0 > w0) ranky++;
          else rankw++;
          if (z0 > w0) rankz++;
          else rankw++;
          var i1, j1, k1, l1; // The integer offsets for the second simplex corner
          var i2, j2, k2, l2; // The integer offsets for the third simplex corner
          var i3, j3, k3, l3; // The integer offsets for the fourth simplex corner
          // simplex[c] is a 4-vector with the numbers 0, 1, 2 and 3 in some order.
          // Many values of c will never occur, since e.g. x>y>z>w makes x<z, y<w and x<w
          // impossible. Only the 24 indices which have non-zero entries make any sense.
          // We use a thresholding to set the coordinates in turn from the largest magnitude.
          // Rank 3 denotes the largest coordinate.
          i1 = rankx >= 3 ? 1 : 0;
          j1 = ranky >= 3 ? 1 : 0;
          k1 = rankz >= 3 ? 1 : 0;
          l1 = rankw >= 3 ? 1 : 0;
          // Rank 2 denotes the second largest coordinate.
          i2 = rankx >= 2 ? 1 : 0;
          j2 = ranky >= 2 ? 1 : 0;
          k2 = rankz >= 2 ? 1 : 0;
          l2 = rankw >= 2 ? 1 : 0;
          // Rank 1 denotes the second smallest coordinate.
          i3 = rankx >= 1 ? 1 : 0;
          j3 = ranky >= 1 ? 1 : 0;
          k3 = rankz >= 1 ? 1 : 0;
          l3 = rankw >= 1 ? 1 : 0;
          // The fifth corner has all coordinate offsets = 1, so no need to compute that.
          var x1 = x0 - i1 + G4; // Offsets for second corner in (x,y,z,w) coords
          var y1 = y0 - j1 + G4;
          var z1 = z0 - k1 + G4;
          var w1 = w0 - l1 + G4;
          var x2 = x0 - i2 + 2.0 * G4; // Offsets for third corner in (x,y,z,w) coords
          var y2 = y0 - j2 + 2.0 * G4;
          var z2 = z0 - k2 + 2.0 * G4;
          var w2 = w0 - l2 + 2.0 * G4;
          var x3 = x0 - i3 + 3.0 * G4; // Offsets for fourth corner in (x,y,z,w) coords
          var y3 = y0 - j3 + 3.0 * G4;
          var z3 = z0 - k3 + 3.0 * G4;
          var w3 = w0 - l3 + 3.0 * G4;
          var x4 = x0 - 1.0 + 4.0 * G4; // Offsets for last corner in (x,y,z,w) coords
          var y4 = y0 - 1.0 + 4.0 * G4;
          var z4 = z0 - 1.0 + 4.0 * G4;
          var w4 = w0 - 1.0 + 4.0 * G4;
          // Work out the hashed gradient indices of the five simplex corners
          var ii = i & 255;
          var jj = j & 255;
          var kk = k & 255;
          var ll = l & 255;
          // Calculate the contribution from the five corners
          var t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0;
          if (t0 < 0) n0 = 0.0;
          else {
            var gi0 = (perm[ii + perm[jj + perm[kk + perm[ll]]]] % 32) * 4;
            t0 *= t0;
            n0 = t0 * t0 * (grad4[gi0] * x0 + grad4[gi0 + 1] * y0 + grad4[gi0 + 2] * z0 + grad4[gi0 + 3] * w0);
          }
          var t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1;
          if (t1 < 0) n1 = 0.0;
          else {
            var gi1 = (perm[ii + i1 + perm[jj + j1 + perm[kk + k1 + perm[ll + l1]]]] % 32) * 4;
            t1 *= t1;
            n1 = t1 * t1 * (grad4[gi1] * x1 + grad4[gi1 + 1] * y1 + grad4[gi1 + 2] * z1 + grad4[gi1 + 3] * w1);
          }
          var t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2;
          if (t2 < 0) n2 = 0.0;
          else {
            var gi2 = (perm[ii + i2 + perm[jj + j2 + perm[kk + k2 + perm[ll + l2]]]] % 32) * 4;
            t2 *= t2;
            n2 = t2 * t2 * (grad4[gi2] * x2 + grad4[gi2 + 1] * y2 + grad4[gi2 + 2] * z2 + grad4[gi2 + 3] * w2);
          }
          var t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3;
          if (t3 < 0) n3 = 0.0;
          else {
            var gi3 = (perm[ii + i3 + perm[jj + j3 + perm[kk + k3 + perm[ll + l3]]]] % 32) * 4;
            t3 *= t3;
            n3 = t3 * t3 * (grad4[gi3] * x3 + grad4[gi3 + 1] * y3 + grad4[gi3 + 2] * z3 + grad4[gi3 + 3] * w3);
          }
          var t4 = 0.6 - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4;
          if (t4 < 0) n4 = 0.0;
          else {
            var gi4 = (perm[ii + 1 + perm[jj + 1 + perm[kk + 1 + perm[ll + 1]]]] % 32) * 4;
            t4 *= t4;
            n4 = t4 * t4 * (grad4[gi4] * x4 + grad4[gi4 + 1] * y4 + grad4[gi4 + 2] * z4 + grad4[gi4 + 3] * w4);
          }
          // Sum up and scale the result to cover the range [-1,1]
          return 27.0 * (n0 + n1 + n2 + n3 + n4);
        }
      };

      function buildPermutationTable(random) {
        var i;
        var p = new Uint8Array(256);
        for (i = 0; i < 256; i++) {
          p[i] = i;
        }
        for (i = 0; i < 255; i++) {
          var r = i + ~~(random() * (256 - i));
          var aux = p[i];
          p[i] = p[r];
          p[r] = aux;
        }
        return p;
      }
      SimplexNoise._buildPermutationTable = buildPermutationTable;

      function alea() {
        // Johannes Baagøe <baagoe@baagoe.com>, 2010
        var s0 = 0;
        var s1 = 0;
        var s2 = 0;
        var c = 1;

        var mash = masher();
        s0 = mash(' ');
        s1 = mash(' ');
        s2 = mash(' ');

        for (var i = 0; i < arguments.length; i++) {
          s0 -= mash(arguments[i]);
          if (s0 < 0) {
            s0 += 1;
          }
          s1 -= mash(arguments[i]);
          if (s1 < 0) {
            s1 += 1;
          }
          s2 -= mash(arguments[i]);
          if (s2 < 0) {
            s2 += 1;
          }
        }
        mash = null;
        return function() {
          var t = 2091639 * s0 + c * 2.3283064365386963e-10; // 2^-32
          s0 = s1;
          s1 = s2;
          return s2 = t - (c = t | 0);
        };
      }
      function masher() {
        var n = 0xefc8249d;
        return function(data) {
          data = data.toString();
          for (var i = 0; i < data.length; i++) {
            n += data.charCodeAt(i);
            var h = 0.02519603282416938 * n;
            n = h >>> 0;
            h -= n;
            h *= n;
            n = h >>> 0;
            h -= n;
            n += h * 0x100000000; // 2^32
          }
          return (n >>> 0) * 2.3283064365386963e-10; // 2^-32
        };
      }
      // common js
      exports.SimplexNoise = SimplexNoise;
      // nodejs
      {
        module.exports = SimplexNoise;
      }

    })();
    });
    var simplexNoise_1 = simplexNoise.SimplexNoise;

    const simplex = new simplexNoise();

    function sum_octave(num_iterations, x, y, scale, persistence, sigmoid_intensity) {
      let noise = 0;
      let maxAmp = 0;
      let amp = 1;
      let freq = 1 / scale;

      for (let i = 0; i < num_iterations; i++) {
        noise += simplex.noise3D(x * freq, y * freq, i) * amp;
        maxAmp += amp;
        amp *= persistence;
        freq *= 2;
      }
      var output = apply_sigmoid(noise / maxAmp, sigmoid_intensity);
      return output;
    }

    function apply_sigmoid(value, intensity) {
      if (intensity === 0) return value;
      return 2 * sigmoid(value * intensity) - 1;
    }

    function sigmoid(x) {
      return 1 / (1 + Math.exp(-x));
    }

    let sketch = function (p) {
      let THE_SEED;
      let noise_grid;
      let slope_grid;
      let pds;
      let points;

      const darkness = 3;

      const scale = 200;
      const persistence = 0.35;
      const sigmoid_intensity = 0;

      const grid_dim_x = 900;
      const grid_dim_y = 900;
      const padding = 40;
      const canvas_dim_x = grid_dim_x + 2 * padding;
      const canvas_dim_y = grid_dim_y + 2 * padding;
      const cell_dim = 2;
      const nx = grid_dim_x / cell_dim;
      const ny = grid_dim_y / cell_dim;

      p.setup = function () {
        p.createCanvas(canvas_dim_x, canvas_dim_y);
        p.pixelDensity(4);
        p.noStroke();
        p.noLoop();
        p.fill(0);

        THE_SEED = p.floor(p.random(9999999));
        p.randomSeed(THE_SEED);
      };

      p.draw = function () {
        p.background('#fd0');
        reset();
        //display();
        displaySample();
      };

      p.keyPressed = function () {
        if (p.keyCode === 80) p.saveCanvas('sketch_' + THE_SEED, 'jpeg');
      };

      function displaySample() {
        //p.fill('#f00');
        p.stroke(0);
        p.strokeWeight(1.5);
        p.translate(padding, padding);
        for (const pnt of points) {
          p.point(pnt[0], pnt[1]);
        }
      }

      function reset() {
        noise_grid = build_noise_grid(0.5, 0);
        slope_grid = build_slope_grid();

        console.log(slope_grid);

        //console.log(slope_grid);
        pds = createPoisson(slope_grid);
        points = pds.fill();
        //console.log(points);

        p.randomSeed(THE_SEED);
      }

      function build_noise_grid(baseline, offset_mag) {
        return [...Array(ny + 1)].map((_, y) =>
          [...Array(nx + 1)].map(
            (_, x) =>
              sum_octave(16, x, y, scale, persistence, sigmoid_intensity) +
              offset_mag * (center_offset(x, y) - baseline)
          )
        );
      }

      function build_slope_grid() {
        const grid = [];
        for (let i = 0; i < noise_grid.length; i++) {
          const row = [];
          for (let j = 0; j < noise_grid[i].length; j++) {
            row.push(map_to_diff(j, i));
          }
          grid.push(row);
        }
        return grid;
      }

      function createPoisson(slope) {
        return new poissonDiskSampling({
          shape: [grid_dim_x, grid_dim_y],
          minDistance: 1,
          maxDistance: 15,
          distanceFunction: function (p) {
            //console.log(p);
            return Math.pow(slope[Math.floor(p[1] / cell_dim)][Math.floor(p[0] / cell_dim)], 8); // value between 0 and 1
          },
        });
      }

      // Output range [-1, 1]
      function center_offset(x, y) {
        return 1 - distance_from_centre(x, y) * 3;
      }

      // Output range: [0, 1] (within tangent circle);
      function distance_from_centre(x, y) {
        return Math.sqrt(Math.pow(nx / 2 - x, 2) + Math.pow(ny / 2 - y, 2)) / nx;
      }

      function map_to_diff(x, y) {
        const slope = get_slope(x, y);
        const mag = Math.sqrt(Math.pow(slope[0], 2) + Math.pow(slope[1], 2));
        const dir = get_dir(...slope);
        const dirdiff = Math.abs(Math.PI - dir);
        return p.constrain(1 - mag * dirdiff * darkness, 0, 1);
      }

      function get_slope(x, y) {
        if (y <= 0 || y >= noise_grid.length - 1) return [0, 0];
        if (x <= 0 || x >= noise_grid[y].length - 1) return [0, 0];
        const n = noise_grid[y - 1][x];
        const s = noise_grid[y + 1][x];
        const w = noise_grid[y][x - 1];
        const e = noise_grid[y][x + 1];

        const lat = s - n;
        const lon = e - w;

        return [lon, lat];
      }

      function get_dir(w, h) {
        let v = p.createVector(w, h);
        return v.heading();
      }
    };
    new p5(sketch);

}));
