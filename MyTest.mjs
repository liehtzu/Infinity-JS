function a(x)
{
    var r = { v1: 0.0, v2: 0.0};
    r.v1 = 2.0 * x;
    r.v2 = x * x;
    return r;
}

function b(x, r)
{
    let changedObj = a(x);
    r.v1 = changedObj.v1;
    r.v2 = changedObj.v2;
    
}

var rr = { v1: 100.0, v2: 1000.0 };
b(5.0, rr);

console.log(rr);
