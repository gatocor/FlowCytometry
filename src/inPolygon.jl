function onSegment(p::Tuple{<:Real,<:Real}, q::Tuple{<:Real,<:Real}, r::Tuple{<:Real,<:Real})
     
    if ((q[1] <= max(p[1], r[1])) &&
        (q[1] >= min(p[1], r[1])) &&
        (q[2] <= max(p[2], r[2])) &&
        (q[2] >= min(p[2], r[2])))
        return true
    else
        return false
    end
end

function orientation(p::Tuple{<:Real,<:Real}, q::Tuple{<:Real,<:Real}, r::Tuple{<:Real,<:Real})
     
    val = (((q[2] - p[2]) *
            (r[1] - q[1])) -
           ((q[1] - p[1]) *
            (r[2] - q[2])))
            
    if val == 0
        return 0
    elseif val > 0
        return 1 # Collinear
    else
        return 2 # Clock or counterclock
    end
end

function doIntersect(p1, q1, p2, q2)
     
    # Find the four orientations needed for 
    # general and special cases
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)
 
    # General case
    if (o1 != o2) && (o3 != o4)
        return true
    end
     
    # Special Cases
    # p1, q1 and p2 are collinear and
    # p2 lies on segment p1q1
    if (o1 == 0) && (onSegment(p1, p2, q1))
        return true
    end
    # p1, q1 and p2 are collinear and
    # q2 lies on segment p1q1
    if (o2 == 0) && (onSegment(p1, q2, q1))
        return true
    end
    # p2, q2 and p1 are collinear and
    # p1 lies on segment p2q2
    if (o3 == 0) && (onSegment(p2, p1, q2))
        return true
    end
    # p2, q2 and q1 are collinear and
    # q1 lies on segment p2q2
    if (o4 == 0) && (onSegment(p2, q1, q2))
        return true
    end

    return false
end

"""
    function isInsideGate(p::Tuple{<:Real,<:Real}, gate::FlowCytometryGate)

Function that given a point and a gate, checks which entries are inside the gate polygon.

**Arguments**:
 - **p::Tuple{<:Real,<:Real}** Point to check if inside polygon.
 - **gate::FlowCytometryGate** Gate to check if the rows are inside the polygon of the gate.

**Returns**:
 Bool. True entries are the ones inside the gate polygon.
"""
function isInsideGate(p::Tuple{<:Real,<:Real}, gate::FlowCytometryGate)
     
    points = gate.polygon
    n = length(points)
     
    # There must be at least 3 vertices
    # in polygon
    if n < 3
        return false
    end
         
    # Create a point for line segment
    # from p to infinite
    extreme = (10E15, p[2])
     
    # To count number of points in polygon
      # whose y-coordinate is equal to
      # y-coordinate of the point
    decrease = 0
    count = 0
    i = 1
     
    while true
        next = i % n + 1
         
        #Avoid lines touching extremal points
        if(points[i][2] == p[2])
            decrease += 1
        end
        if(points[next][2] == p[2])
            decrease += 1
        end
         
        # Check if the line segment from 'p' to 
        # 'extreme' intersects with the line 
        # segment from 'polygon[i]' to 'polygon[next]'
        if (doIntersect(points[i],
                        points[next],
                        p, extreme))
                             
            # If the point 'p' is collinear with line 
            # segment 'i-next', then check if it lies 
            # on segment. If it lies, return true, otherwise false
            if orientation(points[i], p,
                           points[next]) == 0
                return onSegment(points[i], p,
                                 points[next])
            end
                                  
            count += 1
        end
             
        i = next
         
        if (i == 1)
            break
        end
    end
             
    # Reduce the count by decrease amount
      # as these points would have been added twice
    count -= decrease
     
    # Return true if count is odd, false otherwise
    return (count % 2 == 1)
end


"""
    function isInsideGate(x::Matrix{<:Real}, gate::FlowCytometryGate)

Function that given the entries of two channels and a gate checks which entries are inside the gate polygon.

**Arguments**:
 - **x::Matrix{<:Real}** Matrix with N cells and 2 columns (channels).
 - **gate::FlowCytometryGate** Gate to check if the rows are inside the polygon of the gate.

**Returns**:
 Vector{Bool} with N entries. True entries are the ones inside the gate polygon.
"""
function FlowCytometry.isInsideGate(x::Matrix{<:Real}, gate::FlowCytometryGate)

    s = size(x)

    if s[2] != 2
        error("Only matrices with two columns corresponding to two channels are possible.")
    end

    inGate = fill(false,s[1])
    for i in 1:s[1]
        inGate[i] = isInsideGate((x[i,1],x[i,2]),gate)
    end

    return inGate
end