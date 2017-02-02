<template>
  <div class="home">
    <table>
      <tr>
        <td rowspan="2"><img src="static/logo.png"></img>
        </td>
        <td colspan="2" class="home-nsp-intro">Nucleic acids<br>Structure Prediction</td>
      </tr>
      <tr>
        <td><a href="#/doc">Documentation</a></td>
        <td><a href="https://github.com/hust220/nsp">Github</a></td>
      </tr>
    </table>
  </div>
</template>

<script>
  import { bus } from '../bus.js'

  export default {
  }

</script>

<style>
  .home {
    height: 100%;
    padding: 0px;
  }

  .home-nsp-intro {
    font-size: 1.5em;
  }

  .home td {
    padding: 10px 15px;
  }

  .home table {
    margin: 200px auto;
  }

</style>