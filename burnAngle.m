function alpha = burnAngle (dV, V)
    alpha = acosd(dot(dV,V)/(norm(dV)*norm(V)));

    
    
end