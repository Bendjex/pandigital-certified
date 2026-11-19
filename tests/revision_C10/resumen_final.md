# Revisión adversaria C10 — resumen final

(Respuesta del revisor al autor tras las dos comprobaciones finales. El detalle de cada problema está en `informe_revision_C10.md`.)

Las dos comprobaciones que aprobaste terminaron bien y ya están en el informe, que te acabo de reenviar: ninguna destapa un problema nuevo.

**1. Contador de 4.1 en n = 9, cualquier orden, tablas de 5 dígitos.** 947,1 millones de llamadas a `finish_real` y 0 casos con cota superior +∞. El defecto de `code_for`/`finish_real` queda como latente: no afecta a ninguna casilla. Sigue haciendo falta el guard `if (!isfinite(hi)) return D_SING_W;` y tratar ese caso en la prueba del lema 9.

**2. Motor X para δ, n = 8, cualquier orden, con la ventana sin menos unario (1.1).**
- Resultado: 127 separaciones de 127, 492 hits, 0 sin decidir frente a las dos ventanas, 0 valores dentro de la ventana sin menos unario y 0 operaciones de motor B sin decidir, en 20 min de CPU.
- Entradas: el radio sale del óptimo sin menos de `c11.jsonl` (1,2001163304·10⁻¹¹), y el fichero de resolución lo reconstruí aquí con el mismo hash que el publicado.
- Contraste con tus datos: los hits por separación coinciden, las 127, con tus registros publicados de radio estrecho.

La afirmación del artículo resulta por tanto cierta, pero **tu paquete sigue sin respaldarla**. Te falta repetir la ejecución en tu máquina, añadir sus registros a `registros_c5_sesion/n8/`, regenerar `estado_c5.csv` (hoy da esa fila por pendiente aunque `c11.jsonl` ya trae el óptimo) y añadir la entrada de δ a `radii.txt`. Los pasos exactos, mis registros y las trazas están en `pruebas_revisor/delta_n8_sin_menos/` (empieza por `pasos.txt`).

**Veredicto actualizado:**
- Nada de lo encontrado invalida el teorema 1 ni el teorema 2 tal como están enunciados (con sus marcas TODO); los defectos latentes de 1.4 y 4.1 dan 0 casos en los recuentos completos.
- La última frase del teorema 2 para δ, n = 8, cualquier orden es cierta según mi ejecución, pero falta en los datos entregados.

Sigue sin hacer:
- La tercera comprobación opcional (cobertura estructural con conjuntos de hasta 9 dígitos) no la he lanzado porque no la aprobaste; dime si la quieres.
- Las páginas de Friedman y MathWorld no son accesibles desde esta sesión, así que la columna «Previous» de la tabla 1, la tabla 7 y §3.3 quedan sin contrastar.
