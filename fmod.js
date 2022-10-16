Math.fmod = function (a,b)
{
    let x = a / b;
    let y = Math.floor( a / b );
    let z = Math.floor( a / b ) * b;
    let w = a - Math.floor( a / b ) * b;
    console.log( x );
    console.log( y );
    console.log( z );
    console.log( w );
    return Number((a - (Math.floor(a / b) * b)).toPrecision(15));
};

function findMod(a, b)
{
    let mod;
    // Handling negative values
    if (a < 0)
        mod = -a;
    else
        mod =  a;
    if (b < 0)
        b = -b;
 
    // Finding mod by
    // repeated subtraction
     
    while (mod >= b)
        mod = mod - b;
 
    // Sign of result typically
    // depends on sign of a.
    if (a < 0)
        return -mod;
 
    return mod;
}
 

let s = -1.0e-10;

let fm2 = Math.fmod(s, 2.0 );
let fm4 = Math.fmod(s, 4.0 );


let fm22 = findMod(s, 2.0 );
let fm24 = findMod(s, 4.0 );

console.log( fm22 );
console.log( fm24 );
