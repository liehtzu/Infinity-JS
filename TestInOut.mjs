// Function to process in out parameters

import { InOutF } from "./InOut.mjs";

console.log( "Test In Out Parameters" );

var d1 = 1.0;
var d2 = 2.0;
var o1 = { Double: 10.0 };
var o2 = { Double: 20.0 };

console.log( "MAIN: Variables BEFORE function call" );
console.log( "d1 = ", d1 );
console.log( "d2 = ", d2 );
console.log( "o1 = ", o1 );
console.log( "o2 = ", o2 );

// Calling the function
InOutF( d1, d2, o1, o2 );

console.log( "MAIN: Variables AFTER function call" );
console.log( "d1 = ", d1 );
console.log( "d2 = ", d2 );
console.log( "o1 = ", o1 );
console.log( "o2 = ", o2 );
