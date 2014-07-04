<?
session_start();
$sid=session_id();
////////////////////////////////////////////////////////////
//ENTRADAS
////////////////////////////////////////////////////////////
foreach(array_keys($_GET) as $field){
    $$field=$_GET[$field];
}
foreach(array_keys($_POST) as $field){
    $$field=$_POST[$field];
}
$dut=-5;

function sign($n) {
  return $n>=0?+1:-1;
}

function value2dec($ang)
{
  if(preg_match("/:/",$ang)){
    $parts=preg_split("/:/",$ang);
    $ang_deg=$parts[0];
    $ang_min=$parts[1];
    $ang_sec=$parts[2];
  }else if(preg_match("/\s+/",$ang)){
    $parts=preg_split("/\s+/",$ang);
    echo "Espacio.";
    $ang_deg=$parts[0];
    $ang_min=$parts[1];
    $ang_sec=$parts[2];
  }else{
    $ang_deg=$ang;
    $ang_min=0.0;
    $ang_sec=0.0;
  }
  $sgn=sign($ang);
  $abs_ang_dec=abs($ang_deg);
  $ang_dec=$sgn*($abs_ang_dec+
		 $ang_min/60.0+$ang_sec/3600.0);
  return $ang_dec;
}
?>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
</head>
<center>
<img src="images/aristarco.jpg" height=200>
<H1>Ocultación de Marte por la Luna<br/>
Julio 5 de 2014<br/>
</H1>
<H2>
Campaña de Observación: 
<a href=http://bit.ly/aristarco-0705>http://bit.ly/aristarco-0705</a><br/>
Instrucciones Detalladas: 
<a href=http://bit.ly/aristarco-0705-detalles>http://bit.ly/aristarco-0705-detalles</a>
</H2>
</center>

<p>En esta página podrás calcular las condiciones de la ocultación de
Marte por la Luna del próximo 5 de Julio de 2014 para tu ubicación
específica.  Los observadores que se encuentran ubicados sobre la
franja de totalidad (ver mapa abajo) pueden determinar aquí los
tiempos de los contactos teniendo en cuenta o no el efecto de adelanto
del fenómeno debido a la velocidad finita de la luz.</p>

<H3><a name="Tiempos">Tiempos de Contacto</a></H3>

<p>A continuación ingrese su posición geográfica exacta (tal y como es
provista por un GPS).</p>
<form action="?#Tiempos">
Longitud: <input type="text" name="lon" value="<?echo $lon?>"><br/>
<i style=font-size:12px>Formato: +GG.gggg, +GG:MM:SS.sss, +GG MM SS.sss.  Ejemplos: -75 35 23.99, -75:35:23.99, -75.43555</i>
<p></p>
Latitud: <input type="text" name="lat" value="<?echo $lat?>"><br/>
<i style=font-size:12px>Formato: +GG.gggg, +GG:MM:SS.sss, +GG MM SS.sss.  Ejemplos: -75 35 23.99, -75:35:23.99, -75.43555</i>
<p></p>
Uso horario: <input type="text" name="dut" value="<?echo $dut?>"><br/>
<i style=font-size:12px>Formato: +UTC.  Ejemplo: -5, -4.5</i>
<p></p>
Altura: <input type="text" name="alt" value="<?echo $alt?>"> metros<br/>
<i style=font-size:12px>Formato: +AAAA.aaa.  No use "," ni separadores de miles.</i>
<p></p>
<input type='submit' name='accion' value='Calcule'>
</form>
<?
////////////////////////////////////////
//CALCULO
////////////////////////////////////////
if($accion=="Calcule")
{
  //PREPARE INPUTS
  $lon=value2dec($lon);
  $lat=value2dec($lat);

  //CALL PROGRAM
  $out=shell_exec("./contact-times.out $lon $lat $alt $dut 2> tmp/contact-times.log");
  if(!preg_match("/\w/",$out)){
    echo "<i style=color:red>No hay ocultación en esta locación</i>";
  }else{
  $times=preg_split("/,/",$out);
  
echo<<<CONTENIDO
<H4>Observador</H4>
  <b>Longitud:</b>$lon<br/>
  <b>Latitud:</b>$lat<br/>
  <b>Altura:</b>$alt metros<br/>
  <b>Huso horario:</b>UTC$dut<br/>
  <p></p>
  <table border=1>
  <tr>
  <td><b>Contacto</b></td>
  <td><b>Tiempo Estimado</b><br/>(Velocidad de la luz Infinita)</td>
  <td><b>Tiempo Estimado</b><br/>(Velocidad de la luz finita)</td>
  </tr>
  <tr>
  <td>Contacto 1 (Marte toca la Luna)</td>
  <td>$times[7]</td>
  <td>$times[0]</td>
  </tr>
  <tr>
  <td>Contacto 1 punto medio (La mitad de Marte esta ocultada)</td>
  <td>$times[8]</td>
  <td>$times[1]</td>
  </tr>
  <tr>
  <td>Contacto 2 (Marte esta ocultado completamente)</td>
  <td>$times[9]</td>
  <td>$times[2]</td>
  </tr>
  <tr>
  <td>Contacto 3 (Marte empieza a emerger)</td>
  <td>$times[10]</td>
  <td>$times[3]</td>
  </tr>
  <tr>
  <td>Contacto 3 punto medio (La mitad de Marte ha emergido)</td>
  <td>$times[11]</td>
  <td>$times[4]</td>
  </tr>
  <tr>
  <td>Contacto 4 (Fin de la ocultación)</td>
  <td>$times[12]</td>
  <td>$times[5]</td>
  </tr>
  </table>
  <p></p>
  <b>Duración del evento:</b>$times[6] horas
  <p></p>
CONTENIDO;
  //GENERATE CORD
  $filename=sprintf("data/cord-%+.4lf_%+.4lf_%+.3lf.dat",$lon,$lat,$alt);
  shell_exec("MPLCONFIGDIR=/tmp python plot-cords.py $sid $filename &> tmp/cord.log");
  //SHOW CORD
  sleep(1);
echo<<<CONTENIDO
  <b>Trayectoria de Marte Respecto a la Luna:</b>
<p></p>
<img src="tmp/cords-$sid.png" height=400>
CONTENIDO;
  }
}
?>
<H3>Franja de Totalidad</H3>
<center><img src="images/totalidad.png"></center>
<H4><a name="precalculados">Contactos Precalculados</a><br>
Campaña de Observación: 
<a href=http://bit.ly/aristarco-0705>http://bit.ly/aristarco-0705</a>
</H4>

  <center><b>Medellín, Colomia</b>
<center><img src="images/ContactosMedellin.jpg"></center>
<p></p>
  <center><b>Bogotá, Colomia</b>
<center><img src="images/ContactosBogota.jpg"></center>
  <center><b>Cali, Colomia</b>
<center><img src="images/ContactosCali.jpg"></center>
  <center><b>Pasto, Colomia</b>
<center><img src="images/ContactosPasto.jpg"></center>
  <center><b>Desierto de la Tatacoa, Colomia</b>
<center><img src="images/ContactosTatacoa.jpg"></center>

</html>
